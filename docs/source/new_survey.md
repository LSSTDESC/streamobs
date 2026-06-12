# Adding a new photometric survey

To add a new survey:

1. Create a new subdirectory in `data/surveys/` + its release:
   ```bash
   mkdir -p data/surveys/new_survey_name_release
   ```

2. Add magnitude limit maps for each band:
   - Files should be HEALPix format (`.hsp`) or `.fit`
   - Example: `des_y6_g_band_nside_128.hsp`
   - Optional: one may add new completeness and photometric errors data. If not provided, the ones in data/others will be used.

3. Create a corresponding survey configuration in `config/surveys/`:
   ```yaml
   name: new_survey
   release: its_release
   survey_files: 
      # Path to data files (leave empty for default location)
      file_path: ''

      # Band-specific magnitude limit maps
      maglim_map_g: new_survey_maglim_g_band.hsp
      maglim_map_r: new_survey_maglim_r_band.hsp

      # Band-independent maps. Keep by defaults files for completeness, ebv map, and photometric errors
      ebv_map: ebv_sfd98_fullres_nside_4096_ring_equatorial.fits
      completeness: stellar_efficiency_cutr.csv
      log_photo_error: photoerror_r.csv
   ```

4. Update the data archive and upload to Zenodo

5. Document the new survey in this document.