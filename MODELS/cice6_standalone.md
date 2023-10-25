# Table of contents
- [Table of contents](#table-of-contents)
- [Table of runs ](#table-of-runs-)
- [Results ](#results-)
  
# Table of runs <a name="table_of_runs"></a>
| run date  | duration    | key parameters  | ocean forcing | atmosphere forcing | link to ice_in                                                                                         |
| --------- | ----------- | --------------- | ------------- | ------------------ | ------------------------------------------------------------------------------------------------------ |
| 03 Jul 23 | 1yr         | $MLD = 20m$     | AOM2          | JRA55do            |                                                                                                        |
|           | (2005)      | $NDTE = 240$    | (monthly)     | (3hrly)            |                                                                                                        |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 06 Jul 23 | 1yr         | $MLD = nil$     | AOM2          | JRA55do            |                                                                                                        |
|           | (2005)      | $NDTE = 240$    | (monthly)     | (3hrly)            |                                                                                                        |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 08 Jul 23 | 1yr         | $MLD = 50m$     | BRAN          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230708_cice6p4_1yr_bran_mld50m/ice_in                  |
|           | (2005)      | $NDTE = 370$    | (daily)       | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230708_cice6p4_1yr_bran_mld50m/ice_diag.d              |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 08 Jul 23 | 1yr         | $MLD = 50m$     | AOM2          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230708_cice6p4_1yr_ndte3710_mld50m/ice_in              |
|           | (2005)      | $NDTE = 3710$   | (monthly)     | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230708_cice6p4_1yr_ndte3710_mld50m/ice_diag.d          |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 11 Jul 23 | 1yr         | $MLD = 50m$     | BRAN          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_bran_ndte3710_mld50m/ice_in         |
|           | (2005)      | $NDTE = 3710$   | (daily)       | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_bran_ndte3710_mld50m/ice_diag.d     |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 11 Jul 23 | 1yr         | $MLD = 50m$     | BRAN          | ERA5               | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_ndte370_mld50m/ice_in               |
|           | (2005)      | $NDTE = 370$    | (daily)       | (hourly)           | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_ndte370_mld50m/ice_diag.d           |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 11 Jul 23 | 1yr         | $MLD = 50m$     | AOM2          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_era5_bran_mld50m/ice_in             |
|           | (2005)      | $NDTE = 3710$   | (monthly)     | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230711_cice6p4_1yr_era5_bran_mld50m/ice_diag.d         |
|           |             | $trestore = 90$ |               |                    |                                                                                                        |
| 13 Jul 23 | 1yr         | $MLD = 50m$     | AOM2          | ERA5               | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_era5_dailyOC_trestore30/ice_in      |
|           | (2005)      | $NDTE = 370$    | (daily)       | (hourly)           | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_era5_dailyOC_trestore30/ice_diag.d  |
|           |             | $trestore = 30$ |               |                    |                                                                                                        |
| 13 Jul 23 | 1yr         | $MLD = 50m$     | BRAN          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_jra55_bran_trestore30/ice_in        |
|           | (2005)      | $NDTE = 370$    | (daily)       | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_jra55_bran_trestore30/ice_diag.d    |
|           |             | $trestore = 30$ |               |                    |                                                                                                        |
| 13 Jul 23 | 1yr         | $MLD = 50m$     | AOM2          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_jra55_dailyOC_trestore30/ice_in     |
|           | (2005)      | $NDTE = 370$    | (daily)       | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_1yr_jra55_dailyOC_trestore30/ice_diag.d |
|           |             | $trestore = 30$ |               |                    |                                                                                                        |
| 27 Jul 23 | 5yr         | $MLD = 50m$     | AOM2          | JRA55do            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_5yr_baseline/ice_in                     |
|           | (2005-2010) | $NDTE = 240$    | (daily)       | (3hrly)            | https://github.com/dpath2o/AFIM/blob/main/etc/20230713_cice6p4_5yr_baseline/ice_diag.d                 |
|           |             | $trestore = 30$ |               |                    |                                                                                                        |

# Results <a name="results"></a>