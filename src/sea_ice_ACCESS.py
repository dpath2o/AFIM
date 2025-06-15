import os, sys
import intake
import datetime
import xarray as xr
import pandas as pd
import numpy  as np
from pathlib  import Path

class SeaIceACCESS:
    def __init__(self, **kwargs):
        if self.client is None:
            print("Warning: No Dask client was provided to SeaIceACCESS.")
        if self.cice_vars_reqd is None:
            raise ValueError("Missing required variable list: self.cice_vars_reqd")

    def list_access_cice_variables(self, experiment=None):
        """
        Lists all available CICE sea ice variables in the ACCESS-NRI Intake catalog
        for the specified experiment.

        Parameters
        ----------
        experiment : str, optional
            ACCESS experiment name (e.g., "access-om2-01-iaf").
            If not provided, falls back to `self.AOM2_dict['experiment']`.

        Returns
        -------
        List[str]
            Sorted list of unique CICE variable names.
        """
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
        Splits a daily dataset into monthly chunks and writes each to a separate Zarr file.

        INPUTS:
        DS : xarray.Dataset
            Daily model output dataset with a 'time' dimension.
        overwrite : bool
            If True, overwrite existing Zarr directories.

        OUTPUTS:
        Writes: D_out/iceh_YYYY-MM.zarr for each month in DS
        """
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
