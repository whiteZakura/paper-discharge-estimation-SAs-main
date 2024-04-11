
## Code Description

This folder includes scripts for mapping water occurrence frequency from Sentinel-1 (S1) data and calculating water  fraction within SAs from Sentinel-2 (S2) data.

### Scripts Overview:

- `S1WaterOccurrenceFrequency.js`: maps the water occurrence frequency for rivers, saving the results to the GEE Asset.
- `S2WaterFraction.js`: Calculates the water fraction in SAs, with results saved to Google Drive.

## Execution Instructions

To utilize these scripts effectively, follow the guidelines below:

1. Paste the code into the Google Earth Engine Code Editor to execute.
2. When running either the S1 or S2 script, ensure you first paste the corresponding imports at the top of the script.
3. The imports for water occurrence frequency in the S2 script require that you have fully executed the S1 script to obtain necessary data inputs.

