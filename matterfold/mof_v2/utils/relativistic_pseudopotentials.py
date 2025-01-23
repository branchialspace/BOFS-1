# Dalcorso pslibrary
!wget https://people.sissa.it/dalcorso/pslibrary/pslibrary.1.0.0.tar.gz
!tar -xzvf pslibrary.1.0.0.tar.gz
!sudo ln -s /content/bin/ld1.x /bin/ld1.x
!echo "/content/bin/" > /content/pslibrary.1.0.0/QE_path
%cd /content/pslibrary.1.0.0
!chmod +x make_all_ps
!./make_all_ps
