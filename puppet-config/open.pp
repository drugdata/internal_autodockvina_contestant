class internalautodockvinacontestant::open
{
  #epel repo
  exec { 'install_epel':
    command => '/bin/yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm',
    creates => '/etc/yum.repos.d/epel.repo'
  }

  Package { ensure => 'installed' }
  $python_deps    = [ 'python2-pip', 'python-psutil', 'python-virtualenv', 'python-tox', 'pylint', 'python-coverage' ]
  $perl_deps      = [ 'perl-Archive-Tar', 'perl-List-MoreUtils' ]
  $other_packages = [ 'libXft', 'openbabel', 'xorg-x11-xauth', 'screen', 'bzip2', 'which', 'rsync', 'ncftp', 'tk-devel', 'libXScrnSaver' ]
  $pymol_deps     = [ 'subversion', 'gcc', 'gcc-c++', 'kernel-devel', 'python-devel', 'tkinter', 'python-pmw', 'glew-devel', 'freeglut-devel', 'libpng-devel', 'freetype-devel', 'libxml2-devel']
  $mesa_packages  = [ 'mesa-libGL-devel','mesa-libEGL-devel','mesa-libGLES-devel' ]
  $pip_packages   = [ 'argparse','psutil','biopython','xlsxwriter','ftpretty','wheel','flake8','lockfile','easywebdav','d3r', 'numpy' ]
  package { $python_deps: }
  package { $perl_deps: }
  package { $other_packages: }
  package { $mesa_packages: }
  package { $pymol_deps: }

  yumrepo { 'giallu-rdkit':
    baseurl             => 'https://copr-be.cloud.fedoraproject.org/results/giallu/rdkit/epel-$releasever-$basearch/',
    descr               => 'Copr repo for rdkit owned by giallu',
    enabled             => 1,
    gpgcheck            => 1,
    gpgkey              => 'https://copr-be.cloud.fedoraproject.org/results/giallu/rdkit/pubkey.gpg',
    repo_gpgcheck       => 0,
    skip_if_unavailable => 'true'
  }

  package { 'python-rdkit':
    require => Yumrepo['giallu-rdkit']
  }

  package { $pip_packages:
    ensure   => 'installed',
    provider => 'pip',
    require  => Package['python2-pip'],
  }

  # manual INSTALL Chimera
  exec { 'install_chimera':
    path => [ '/usr/bin', '/usr/sbin', '/bin', '/usr/local/bin'],
    command => 'cd /tmp;
                chmod a+x chimera*.bin;
                echo "/opt/chimera" | ./chimera*.bin;
                cd /tmp;
                /bin/rm -rf chimera*.bin',
    onlyif => '/bin/test -e /tmp/chimera*.bin',
    creates => '/opt/chimera/bin/chimera'
  }

   # install autodockvina
   exec { 'install_autodockvina':
     path => [ '/usr/bin', '/usr/sbin', '/bin', '/usr/local/bin'],
     command => 'cd /tmp;
                 wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz;
                 tar -zxf autodock_vina_1_1_2_linux_x86.tgz;
                 mv autodock_vina_1_1_2_linux_x86/bin/* /usr/local/bin/.;
                 cd /tmp;
                 /bin/rm -rf autodock_vina_1_1_2_linux_x86*',
     creates => '/usr/local/bin/vina'
   }

   # install mgl tools
   exec { 'install mgltools':
     path => [ '/usr/bin', '/usr/sbin', '/bin', '/usr/local/bin'],
     command => 'cd /tmp;
                 wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz;
                 tar -zxf mgltools_x86_64Linux2_1.5.6.tar.gz;
                 cd mgltools_x86_64Linux2_1.5.6;
                 ./install.sh -d /usr/local/mgltools;
                 cd /tmp;
                 /bin/rm -rf mgltools_x86_64Linux2_1.5.6*',
     creates => '/usr/local/mgltools'
   }
}
# Run the class
class { 'internalautodockvinacontestant::open': }
