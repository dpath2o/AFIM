from __future__ import annotations
class SeaIceACCESS:
    """
    Interface to the ACCESS-NRI Intake catalog for loading ACCESS-OM* CICE daily output
    and writing it into an AFIM-friendly monthly Zarr layout.

    This class is intended to:
      1) Query the ACCESS-NRI Intake catalog for available CICE variables for a given
         ACCESS experiment.
      2) Load a daily dataset (1day frequency) for a selected set of variables and
         time window.
      3) Standardise basic indexing/coordinates (e.g., ensure `nj`/`ni` coordinates
         exist) so downstream merging/processing is reliable.
      4) Persist the daily time series to monthly Zarr groups (one directory per month).

    Notes
    -----
    - The Intake catalog must be installed and the ACCESS-NRI catalog must be available
      as `intake.cat.access_nri`.
    - This class assumes CICE daily output is labelled as a daily mean associated with
      the *following* day (a common convention); `write_ACCESS_to_monthly_zarr()` shifts
      timestamps back by 1 day to align with the intended averaging period.
    - Dask is optional, but recommended for performance. If `self.client` is not set,
      operations may still work but can be slower and may run locally.

    Expected Attributes
    -------------------
    logger : logging.Logger
        Logger used for status and warnings.
    client : dask.distributed.Client or None
        Dask client (optional). If None, a warning is emitted at initialisation.
    AOM2_dict : dict
        Configuration dictionary expected to contain at least `experiment` (str) when
        `experiment` is not passed explicitly to methods.
    dt0_str, dtN_str : str
        Default ISO date strings (e.g. "1993-01-01") used when method arguments are None.
    cice_vars_reqd : list[str]
        Required CICE variable names you want to extract from the loaded dataset.
    cice_var_list : list[str]
        Variable names used to query the catalog. This can match `cice_vars_reqd`, or be
        a broader list (e.g., including coordinates/auxiliary vars).
    D_zarr : str or pathlib.Path
        Output directory where monthly Zarr groups are written.

    Raises
    ------
    ValueError
        If `self.cice_vars_reqd` is not defined at initialisation.
    """

    def __init__(self, **kwargs):
        """
        Initialise the ACCESS-NRI CICE catalog helper.

        Parameters
        ----------
        **kwargs
            Implementation-dependent configuration parameters. Typical usage is to pass
            or set attributes such as `logger`, `client`, `AOM2_dict`, `dt0_str`, `dtN_str`,
            `cice_vars_reqd`, `cice_var_list`, and `D_zarr`.

        Notes
        -----
        - If no Dask client is provided (`self.client is None`), a warning is printed.
        - `self.cice_vars_reqd` must be set before completion; otherwise a ValueError is raised.

        Raises
        ------
        ValueError
            If `self.cice_vars_reqd` is None.
        """
        if self.client is None:
            print("Warning: No Dask client was provided to SeaIceACCESS.")
        if self.cice_vars_reqd is None:
            raise ValueError("Missing required variable list: self.cice_vars_reqd")

    def list_access_cice_variables(self, experiment=None):
        """
        List available daily CICE variables in the ACCESS-NRI Intake catalog.

        This method queries the ACCESS-NRI catalog for the specified experiment and returns
        the unique set of variable names found at daily frequency (frequency="1day").

        Parameters
        ----------
        experiment : str, optional
            ACCESS experiment identifier (e.g., "access-om2-01-iaf").
            If None, the method falls back to `self.AOM2_dict['experiment']`.

        Returns
        -------
        list[str]
            Sorted list of unique CICE variable names available for the experiment at
            1-day frequency. Returns an empty list if no matching results are found.

        Raises
        ------
        ValueError
            If `experiment` is not provided and `self.AOM2_dict['experiment']` is missing.
        RuntimeError
            If the ACCESS-NRI catalog is not available as `intake.cat.access_nri`.

        Notes
        -----
        - The ACCESS-NRI Intake catalog must be installed and configured (i.e., the catalog
          is discoverable under `intake.cat`).
        - Catalog entries sometimes store variable names as a list per asset. This method
          flattens those lists into a unique set.
        """
        import intake
        experiment = experiment or self.AOM2_dict.get('experiment', None)
        if experiment is None:
            raise ValueError("Experiment name must be provided or set in self.AOM2_dict['experiment'].")

        if not hasattr(intake.cat, "access_nri"):
            raise RuntimeError("ACCESS-NRI catalog not found in Intake. Run `intake catalog add access_nri ...` first.")

        self.logger.info(f"Querying ACCESS-NRI catalog for available CICE variables in '{experiment}'...")
        results = intake.cat.access_nri[experiment].search(frequency="1day")

        if results is None or results.df.empty:
            self.logger.warning(f"No CICE results found for experiment '{experiment}' at 1-day frequency.")
            return []

        # Flatten list of lists to extract individual variable names
        var_lists = results.df["variable"].dropna().tolist()
        flat_vars = {var for sublist in var_lists for var in (sublist if isinstance(sublist, list) else [sublist])}

        var_list = sorted(flat_vars)
        self.logger.info(f"Found {len(var_list)} unique variables.")
        return var_list

    def load_ACCESS_OM_CICE(self, experiment=None, dt0_str=None, dtN_str=None):
        """
        Load daily CICE output from the ACCESS-NRI Intake catalog as an xarray.Dataset.

        The method searches the ACCESS-NRI catalog for a daily (frequency="1day") CICE dataset
        matching `self.cice_var_list`, opens it via `to_dataset_dict()`, extracts the variables
        listed in `self.cice_vars_reqd` (if present), and subsets the time dimension to the
        requested window.

        Parameters
        ----------
        experiment : str, optional
            ACCESS experiment identifier. If None, defaults to `self.AOM2_dict['experiment']`.
        dt0_str : str, optional
            Inclusive start date string (ISO-like), e.g. "1993-01-01".
            If None, defaults to `self.dt0_str`.
        dtN_str : str, optional
            Inclusive end date string (ISO-like), e.g. "1993-12-31".
            If None, defaults to `self.dtN_str`.

        Returns
        -------
        xr.Dataset
            Dataset containing the extracted variables on the time subset. The dataset is
            merged from individual variables to ensure consistent coordinate labelling.

        Raises
        ------
        RuntimeError
            If the ACCESS-NRI catalog is not available as `intake.cat.access_nri`.
        KeyError
            If `experiment` is None and `self.AOM2_dict['experiment']` is missing.
        ValueError
            If the catalog search returns zero or multiple daily datasets (expects exactly one).

        Warns
        -----
        UserWarning (via logger)
            If any of `self.cice_vars_reqd` are missing from the loaded dataset.

        Notes
        -----
        - Uses `use_cftime=True` for robust non-standard calendars.
        - Adds explicit `nj` and `ni` coordinates (integer labels) to each extracted variable
          to avoid merge failures or downstream ambiguity when datasets lack labelled indices.
        - The catalog search uses `self.cice_var_list`; variables actually returned are then
          filtered to `self.cice_vars_reqd` (if present).
        """
        import intake
        import numpy as np
        import xarray as xr
        experiment = experiment if experiment is not None else self.AOM2_dict['experiment']
        dt0_str    = dt0_str    if dt0_str    is not None else self.dt0_str
        dtN_str    = dtN_str    if dtN_str    is not None else self.dtN_str
        t_slice    = slice(dt0_str, dtN_str)
        if not hasattr(intake.cat, "access_nri"):
            raise RuntimeError("ACCESS-NRI catalog not found in Intake. Did you add it with `intake catalog add`?")
        result = intake.cat.access_nri[experiment].search(variable=self.cice_var_list, frequency="1day")
        ds_dict = result.to_dataset_dict(xarray_open_kwargs={"use_cftime": True,
                                                             "decode_coords": True,
                                                             "decode_timedelta": False})
        if len(ds_dict) != 1:
            raise ValueError(f"Expected one daily dataset, got {len(ds_dict)}.\n"
                            f"Available keys: {list(ds_dict.keys())}")
        ds = list(ds_dict.values())[0]        
        missing_vars = [var for var in self.cice_vars_reqd if var not in ds]
        if missing_vars:
            self.logger.warning(f"Missing variables in dataset: {missing_vars}")
        extracted_vars = [var for var in self.cice_vars_reqd if var in ds]
        self.logger.info(f"Extracting variables from ESM datastore: {extracted_vars}")
        DS = {var: ds[var].sel(time=t_slice) for var in extracted_vars}
        #DS = {var: ds[var].sel(time=t_slice) for var in self.cice_vars_reqd if var in ds}
        # Add nj, ni coordinate labels so merge works cleanly
        for var, da in DS.items():
            DS[var] = da.assign_coords(nj=("nj", np.arange(da.shape[1])),
                                       ni=("ni", np.arange(da.shape[2])))
        return xr.merge(DS.values())

    def write_ACCESS_to_monthly_zarr(self, DS, overwrite=False):
        """
        Write a daily ACCESS-CICE Dataset into monthly Zarr groups.

        The input dataset is assumed to have a daily time axis. The method:
          1) Creates/ensures the output directory exists (`self.D_zarr`).
          2) Shifts timestamps back by one day to correct the common "daily mean labelled
             by the following day" convention in some CICE outputs.
          3) Splits the dataset into monthly subsets (YYYY-MM).
          4) Writes each month to `iceh_YYYY-MM.zarr` under `self.D_zarr`.

        Parameters
        ----------
        DS : xr.Dataset
            Daily CICE dataset with a `time` dimension. Expected to also have spatial
            dimensions `nj` and `ni`.
        overwrite : bool, default False
            If True, existing monthly Zarr directories are overwritten. If False, existing
            outputs are skipped.

        Returns
        -------
        None
            This method writes Zarr stores to disk and does not return a value.

        Side Effects
        ------------
        Writes monthly Zarr groups to:
            <self.D_zarr>/iceh_YYYY-MM.zarr

        Raises
        ------
        OSError
            If the output directory cannot be created.
        Exception
            Any exception raised by `xarray.Dataset.to_zarr()` will be logged and propagated
            unless caught (current implementation logs errors per month).

        Notes
        -----
        - Time shift: `time = time - 1 day`. If your upstream data is already correctly
          labelled, remove or disable this adjustment.
        - Chunking: this method applies a default chunking
          `{"time": -1, "nj": 540, "ni": 1440}`. Adjust to match your grid shape and
          storage/performance requirements.
        """
        import datetime
        import numpy as np
        import pandas as pd
        import xarray as xr
        from pathlib import Path
        D_out = Path(self.D_zarr)
        D_out.mkdir(parents=True, exist_ok=True)
        DS["time"] = xr.DataArray([t - datetime.timedelta(days=1) for t in DS["time"].values], dims="time", coords={"time": DS["time"].values})
        self.logger.info(f"First timestamp now: {str(DS.time.values[0])}")
        self.logger.info("Shifted all timestamps back by one day to correct for CICE daily average label.")
        time_vals = pd.to_datetime([t.isoformat() for t in DS.time.values])
        mo_strs = np.unique(time_vals.strftime("%Y-%m"))
        for m_str in mo_strs:
            ds_month = DS.sel(time=time_vals.strftime("%Y-%m") == m_str)
            if ds_month.time.size == 0:
                self.logger.warning(f"No data found for {m_str}")
                continue
            P_zarr = Path(D_out,f"iceh_{m_str}.zarr")
            if P_zarr.exists() and not overwrite:
                self.logger.info(f"Skipping {P_zarr} (already exists)")
                continue
            self.logger.info(f"Writing {m_str} with {ds_month.time.size} days → {P_zarr}")
            ds_month = ds_month.chunk({"time": -1, "nj": 540, "ni": 1440})
            try:
                ds_month.to_zarr(P_zarr, mode="w", consolidated=True)
            except Exception as e:
                self.logger.error(f"Failed to write {P_zarr}: {e}")
