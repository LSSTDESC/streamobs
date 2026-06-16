# Update StreamObs data base

If you need to add or modify data files:


1. **Create a new zip archive**:
   ```bash
   cd /path/to/streamobs
   zip -r data.zip data/ \
       -x "*.DS_Store" \
       -x "*__MACOSX*" \
       -x "data/.git*" \
       -x "data/__pycache__/*" \
       -x "*.backup" \
       -x "*.bak" \
       -x "*~"
   ```
   
2. **Upload to Zenodo**:
   - Go to https://zenodo.org/records/18298544
   - Create a new version
   - Upload the `data.zip` file
   - Add release notes describing changes
   - Publish and note the new record ID

3. **Update the download script**:
   - Edit `bin/download_data.py`
   - Update `BASE_DATA_URL` if record ID changed
   - Update `ARCHIVE_SIZE_MB` if size changed
   - Update version in this document

4. **Commit and push**:
   ```bash
   git add .
   git commit -m "Update data archive to version X.Y"
   git push
   ```
