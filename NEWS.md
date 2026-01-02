# phiper 0.2.0

## Minor changes

- Removed generic and S3 methods for the expand_phip_data
- Renamed the internal helpers involved in the data import to match the naming
convention .ph_<function_name>
- documented the internal helpers involved in the data import stage

## Major changes

- Removed the backend argument entirely from the phip_convert, 
phip_convert_legacy, new_phip_data, .resolve_paths and other, less imporatant 
helpers. DuckDB is now the only supported backend
- Changed the structure of the repo. Now all functions related to the 
standard/legacy import workflows live in the standard-conver.R or
legacy-convert.R
- removed the resolve-paths.R file. Moved the functions to utils.R

## Documentation

- Regenerated the docs.


# phiper 0.1.0

## Minor changes

- None.

## Major changes

- None.

## Features

- None.

## Bug fixes

- None.

## Documentation

- None.

## Internal

- None.
