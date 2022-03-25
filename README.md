### Q_robot environment: 

Add the following four lines to your ~/.bashrc file 

```
### Q_robot 
export ROBOT=$HOME/bin/q-robot
export PATH=$PATH:$ROBOT/actions:$ROBOT/friends/vtst/vtstscripts-937
export PYTHONPATH=$PYTHONPATH:$ROBOT/brain
```

### Dependencies

pdftotext
scholarly
selenium
RDkit 
openbabel
xyz2mol 

### Useful Commands

#### To install pip3 

```bash
sudo apt-get update
sudo apt install python3-pip
```

#### To install pdftotext

 https://github.com/jalan/pdftotext

```bash
sudo apt-get install build-essential libpoppler-cpp-dev pkg-config python-dev
pip3 install pdftotext --user
```

#### To install scholarly, selenium

https://github.com/OrganicIrradiation/scholarly

https://github.com/SeleniumHQ/selenium/ 

```bash
pip3 install scholarly --user
pip3 install selenium --user
```

#### RDkit, Openbabel

See website.

#### xyz2mol

https://github.com/jensengroup/xyz2mol# q-robot
