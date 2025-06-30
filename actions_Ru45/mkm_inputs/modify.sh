sed -i 's/thermo: surface-lateral-interaction/thermo: ideal-surface/g' thermo.yaml
sed -i '31i   adjacent-phases: [gas]' thermo.yaml
sed -i '31s/^/  /g' thermo.yaml
