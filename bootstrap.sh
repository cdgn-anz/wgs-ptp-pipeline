#/usr/bin/env bash

case $(id -u) in
    0) 
         echo first: running as root
         echo doing the root tasks...
         sudo apt update && sudo apt install -y build-essential
         echo resizing partition
          (
          echo n   # Create new partition
          echo p   # Primary partition
          echo 2   # Second partition
          echo     # Default first
          echo     # Default last
          echo t   # Change partition ID
          echo 2   # Select 3rd partition
          echo 8e  # Linux LVM
          echo w   # Write partition table
          ) | sudo fdisk /dev/sda
          sudo pvcreate /dev/sda2
          sudo vgextend vagrant-vg /dev/sda2
          sudo lvextend -l +100%FREE /dev/vagrant-vg/root
          sudo resize2fs /dev/vagrant-vg/root
         echo resizing complete
         sudo -u vagrant -i $0  # script calling itself as the vagrant user
         ;;
    *) 
         # running as vagrant user
         wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
         sh Miniconda3-latest-Linux-x86_64.sh -b -p /home/vagrant/miniconda3
         rm Miniconda3-latest-Linux-x86_64.sh
         echo "export PATH=/home/vagrant/miniconda3/bin:/usr/local/bin:$PATH" >> /home/vagrant/.bashrc
         export PATH=/home/vagrant/miniconda3/bin:$PATH
         pip install cython
         pip install git+https://github.com/pytries/datrie.git
         pip install snakemake
         pip install psutil
         pip install pandas
         echo "starting up..."
         # need to put things in $HOME, because can't mount singularity images
         # in /vagrant --- this is a known issue with no fix.
         mkdir -p $HOME/qc_folder/singularity 
         cp /vagrant/startup.sh $HOME/qc_folder/startup.sh
         cp -r /vagrant/pipelines $HOME/qc_folder/
         cd $HOME/qc_folder && bash startup.sh
         echo "All done... Run: snakemake --use-singularity..."
         ;;
esac