#!/bin/bash
MY_PWD=$(pwd)

cd "NetFVFE/two_vessels_angio"
python3 angio.py run

cd $MY_PWD
cd "NetFCFVFE/two_vessels_angio"
python3 angio.py run

cd $MY_PWD
cd "NetFVFEExp/two_vessels_angio"
python3 angio.py run
