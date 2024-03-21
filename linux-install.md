**Linux Install Extra Steps**
*Modify linux-environment.yml to remove 2 lines:*
`  - dss_python=0.10.7.post1`
`  - opendssdirect.py=0.6.1`
*Setup conda environment*
`make init-linux`

*Install via pip*
`pip install dss_python==0.10.7.post1 opendssdirect.py==0.6.1`

*Install python build essentials*
`apt install build-essential python3-dev`


