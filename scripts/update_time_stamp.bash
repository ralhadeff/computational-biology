# this script will run over all folders in the current path (not recursive) and update all folder's timestamps to the newest file within the folder

# this helps fix the issue where sometimes folder timestamps are being updated to the latest time a file was accessed rather than the latest time a file was modified

for folder in $(ls -d */) ; do file=$(ls -rt $folder | tail -1) ; touch -r $folder/$file $folder ; done
