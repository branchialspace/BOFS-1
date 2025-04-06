# Dalcorso fr jpaw generate
!git clone https://github.com/dalcorso/pslibrary.git
!sed -i "s|PWDIR='/path_to_quantum_espresso/'|PWDIR='/content/qe-7.4.1'|" /content/pslibrary/QE_path
!cd /content/pslibrary/rel-pbe && . ../make_ps

# ONCV fr
!git clone https://github.com/pipidog/ONCVPSP.git
!git clone https://github.com/MarioAndWario/ONCVPseudoPack.git

# Dalcorso fr jpaw broken
!wget https://people.sissa.it/dalcorso/pslibrary/pslibrary.1.0.0.tar.gz
!tar -xzvf pslibrary.1.0.0.tar.gz
!sudo ln -s /content/bin/ld1.x /bin/ld1.x
!echo "/content/bin/" > /content/pslibrary.1.0.0/QE_path
%cd /content/pslibrary.1.0.0
!chmod +x make_all_ps
!./make_all_ps

# SSSP library nr
!gdown 1w--QOWnmlPZDy9qJp0NmFakRKHkY8_pt
!gdown 1FoHw9CJT78LItQuaxvrsmgpJ5nrpK1gp
!mkdir -p pseudo_sssp
!tar -xf SSSP_1.3.0_PBE_efficiency.tar.gz -C pseudo_sssp
