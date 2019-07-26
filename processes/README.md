# Initial processing of the raw fast5 data<br>
### Deepbinner
This step uses Deepbinner prior to basecalling, as found [on the github page](https://github.com/rrwick/Deepbinner)
The main command used here is
```
deepbinner realtime --in_dir fast5s --out_dir demultiplexed_fast5s --native
```
This command should move fast5 files from the in-directory to the out-directory, watching the in-directory for any incoming files and acting in realtime. 
### ont_fast5_api
This step uses the ont_fast5_api software from ONT, as found [on their github page](https://github.com/nanoporetech/ont_fast5_api)
Here, we are converting single-read fast5 files to multi-read fast5 files using the code from
```
single_to_multi_fast5
```

