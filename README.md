# Palisades Tahoe rain-on-snow, Feb 24–26, 2026

How much of the snowpack runoff during the storm was rain, and how much was melt? Raw data, model code and a Colab notebook for SNOTEL 784 (Palisades Tahoe, 2,441 m).

## Raw data

Files in `raw/` are exactly as downloaded (SHA-256 below). All times PST (UTC−8) except NASA POWER (UTC).

| File | Source | Retrieved (UTC) | SHA-256 |
|---|---|---|---|
| `NRCS_SNOTEL784_PalisadesTahoe_hourly_2026-02-19_2026-02-28.csv` | [NRCS Report Generator](https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/hourly/784:CA:SNTL/2026-02-19,2026-02-28/TOBS::value,TOBS::qcFlag,TOBS::qaFlag,PREC::value,PREC::qcFlag,PREC::qaFlag,WTEQ::value,WTEQ::qcFlag,WTEQ::qaFlag,SNWD::value,SNWD::qcFlag,SNWD::qaFlag,WSPDV::value,WSPDV::qcFlag,WSPDV::qaFlag,WDIRV::value,RHUM::value,SRADV::value), SNOTEL 784. Provisional data. | 2026-11-08 16:43 | `bf0707a243e984f15de94a5d8259a31a80a806e334c85afadd756c1bba2cb459` |
| `NASA_POWER_hourly_SW_LW_2026-02-19_2026-02-28_UTC.csv` | [NASA POWER hourly API](https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=ALLSKY_SFC_SW_DWN,ALLSKY_SFC_LW_DWN&community=RE&longitude=-120.26&latitude=39.19&start=20260219&end=20260228&format=CSV&time-standard=UTC), CERES SYN1deg all-sky SW and LW down | 2026-11-08 16:56 | `7d6ec3a3bb3b3f9124a0b2332210c5508219ecaac803f4abb3c870079dfe1122` |
| `WesternWeatherGroup_Belmont_BMT_hourly_2026-02-19_2026-02-28_PST.csv` | [Western Weather Group](https://palisadestahoe.westernweathergroup.com/), station BMT (Belmont, beside SNOTEL 784). Unreviewed provider data. | 2026-11-08 17:11 | `849f8b60e5fab56a74e1ed89ab41f0a4cebf513404af2db3eed52e3972fa6849` |
| `MODIS/MODIS_MOD10A1_pixel_SNOTEL784.csv` | MOD10A1.061 ([doi:10.5067/MODIS/MOD10A1.061](https://doi.org/10.5067/MODIS/MOD10A1.061)), pixel containing the station (tile h08v05, row 194, col 1628) and its 3×3 neighbours, extracted from 13 HDF granules (names in the `file` column) | ordered 2026-09-20 | `a7b1bba7ccc76dc62ebcad591a6e9fdb166d42c610b4cab9512aa443319aa7cf` |
