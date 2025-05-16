# io

## read
- `rds` use qs::qread, readr::read_rds
- `csv`, `tsv` use  data.table::fwrite
- Binary `csv` use fst::read_fst
- `xlsx` use readxl::read_xlsx
- `json` use jsonlite::fromJSON
- `yaml` use yaml::read_yaml
- `feather` use arrow::read_feather

## write
- `rds` use qs::qsave
- `csv`, `tsv` use data.table::fwrite
- Binary `csv` use fst::write_fst
- `xlsx` use  writexl::write_xlsx
- `json` use jsonlite::toJSON
- `yaml` use yaml::write_yaml
- `feather` use arrow::write_feather
