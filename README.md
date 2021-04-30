# jm_skyview
Script to retrieve cutout images from NASA Skyview service from fields identified by source name or coordinates.

Some of the functions are adapted from https://github.com/cosmicpudding/skyviewbot/tree/master/skyviewbot

## Usage
```bash
python jm_skyview.py {source_name}
python jm_skyview.py -s {survey} -p {mypath} {source_name}
python jm_skyview.py -s {survey} -p {mypath} {ra_deg},{dec_deg}
```
