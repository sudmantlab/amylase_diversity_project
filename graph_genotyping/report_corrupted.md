# Report on corrupted files 

## CRAM Grch38 1KGP Folder
- Folder path: `/global/scratch/p2p3/pl1_sudmant/human_diversity/1KG30X/data/ERR*/`
- Total number of individuals: `2,504`
- Number of individuals with different md5 hash: `95`

## Additional 698 Related Folder
- Folder path: `/global/scratch/p2p3/pl1_sudmant/human_diversity/1KG30X/additional_698_related/data/ERR398*`
- Total number of individuals: `698`
- Number of individuals with different md5sum: `31`

## Actions Taken
-The corrupted files were re-downloaded and checked.
- The files are now stored in this folder: `/global/scratch/p2p3/pl1_sudmant/alessandroraveane/re_dwn_1kgp/`
- The files were re-indexed.

## Commands Used

To download using the dtn node:
```
cat <ftp_notmatching>.tsv | xargs -P 20 -I{} sh -c 'wget -r -t 0 -c -np -nH --cut-dirs=3 {}'
```
To symlink in the folder `/global/scratch/p2p3/pl1_sudmant/alessandroraveane/re_dwn_1kgp/`:

```
ln -s /global/scratch/users/alessandroraveane/dwn_missdon/<path_to_folder>/ERR*
```

To index in the interactive node / job:

```
ls ERR*/*cram | xargs -P 7 -I{} samtools index {} -M -o {}.crai
```

Both the CRAM Grch38 1KGP Folder and Additional 698 Related Folder had corrupted files, and the same actions were taken to correct the issues. The corrupted files were re-downloaded and checked, and then stored in a new folder with a different name. The newly downloaded files were re-indexed using the samtools index command. 

