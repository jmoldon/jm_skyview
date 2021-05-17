# jm_skyview
Script to retrieve cutout images from NASA Skyview service from fields identified by source name or coordinates and plot them using APLpy.

Some of the functions are adapted from https://github.com/cosmicpudding/skyviewbot/tree/master/skyviewbot

## Usage
```bash
python jm_skyview.py {source_name}
python jm_skyview.py -s {survey} -p {mypath} {source_name}
python jm_skyview.py -s {survey} -p {mypath} {ra_deg},{dec_deg}
```

## Example
```bash
python jm_skyview.py M31  -s DSS --fov 2 -p example
```
[](./example/M31_DSS_2.0d.png)


