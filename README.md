# Data analysis project on GenomeDK

This repository serve as a walk-though of how to begin a data analysis project on the GenomeDK cluter. Be that a PiB, BSc, Msc, or PhD project. Some basic experience with Python programmingis and executing commands in a "terminal" is assumed.

The cluster is a very largre collection of computers with a shared file system. It does not have a screen and a keyboard you can go use. However, by connecting to the cluster from your own computer, you can create and edit files much like if they were on your own machine. The goal here is to take you thrugh all the steps required to make this possible.

## Install Python on your own machine

If you have not done so already, you shuold a distribution of Python called *Anaconda*. Anaconda is simply an easy way of installing Python on Windows, macOS (Mac), and Linux, but it comes with the conda package management system (see below). To install Anaconda, head to [this](https://www.anaconda.com/download) site. Scroll down a bit and click the big green Download button where it says **Python 3.7 version**. When the download has completed, you should follow the platform specific instructions:

Install Anaconda python your local machine.  If you are on a Windows machine you also need to tick the box to add Anaconda python to your PATH when prompted.

* **For Windows:** Double-click the `.exe` file you just downloaded and follow the instructions on the screen. Say yes when asked if want to install Visual Studio Code. The installer will also ask you if you want to download and install a program called Visual Studio Code. Do that too.
* **For OSX:** Double-click the `.pkg` file you just downloaded and follow the instructions on the screen. Make a default installation. The installer will also ask you if you want to download and install a program called Visual Studio Code. Do that too.

## The Terminal

If you are on a Mac or Limux machine, you can use to default Terminal application. If you are on a Windwos macinne, you need to download and install the newest version of Powershell. You can [download it here](https://github.com/PowerShell/PowerShell/releases/download/v7.0.0/PowerShell-7.0.0-win-x64.msi). Double click the installer and. Now open the newly installed Powershell and run: `conda init powershell`. The adantage of PowerShell it that it works the same way as a Linux or OXS terminal. There are some exelent quick tutorials online that introduces the most basic commands in a Linux terminal. When ever I refer to "the terminal" below, it means PoweShell if you are on windows, and the Terminal app if you are on Mac.

## Conda

You need to install packages and programs for use in your analyses and pipelines. Sometimes, however, the versions of packages you need for one project conflicts with the versions you need for other projects that you work on in parallel. Such conflicts seem like an unsolvable problem. Would it not be fantastic if you could create a small world, insulated from the rest of your Anaconda installation. Then that small world could only contain the packages you needed for a single project. If each project had its own isolated world, then there would be no such conflicts. Fortunately, there is a tool that lets you do just that, and its name is Conda. The small worlds that Conda creates are called "environments," and you can create as many as you like, and then switch between them as you switch between your bioinformatics projects. Conda also downloads and installs the packages for you and makes sure that the packages you install in each environment are compatible.  It even makes sure that packages needed by packages  (dependencies) are installed. Conda lets you install mutually compatible versions of software and libraries in an enviromment for your project. By creating an enviromnet for each project, the libraries installed for each project do not interfere. 

## Create an environment on your local machine

When you install Anaconda or Miniconda, Conda makes a single base environment for you. It is called "base" and this is why it says "(base)" at your terminal prompt. You need a conda enviromnet for you project on both your local machine and on the cluster. Lets call both of them 'bircproject' (you can call it anything you like).

The environmnet on your local machine does not need a lot of packages since it mainly serve to let you connect to the cluster. This creates the enviromnet and installs `openssl` and `slurm-jupyter` from my conda chanel:

    conda create -n bircproject -c kaspermunch slurm-jupyter

Say yes (press Enter) when asked to install packages.

## Getting an account on the cluster

The cluster is called GenomeDK and has its own [website](https://genome.au.dk) with lots of information and ducumentation.To get an account on the cluster, you [request one here](https://genome.au.dk/docs/getting-started/#request-access). Below `username` will represent your user name.

On the cluster you have a home folder that only you have access to. That is where end up when you log in. Collaborative projects or projects that use or generate a lot of data projects belong in project folders. If you do a project, we will set up a dedicated project folder for this. 

## Connecting to the cluster

You connect to the cluster from the terminal by executing this command:

    ssh usernmae@login.genome.au.dk

When you do, you are promted for you password to the cluster. You must set up your ssh connection to the cluster so you can connect securely without typing the password every time. First see if you have these two authentication files on your local machine (you can do so by running `ls -a ~/.ssh` in the terminal):

    ~/.ssh/id_rsa
    ~/.ssh/id_rsa.pub

if not, you generate a pair of authentication keys like this. Do not enter a passphrase when prompted - just press enter:

    ssh-keygen -t rsa

Now use `ssh` to create a directory `~/.ssh` on the cluster (assuming your username on the cluster is 'username'):

    ssh username@login.genome.au.dk mkdir -p .ssh

Finally append the public ssh key on your local machine to the file `.ssh/authorized_keys` on the cluster and enter the password one last time (replace username with your cluster user name):

    cat ~/.ssh/id_rsa.pub | ssh username@login.genome.au.dk 'cat >> .ssh/authorized_keys'

From now on you can log into the cluster from your local machine without being prompted for a password.

## Install Python on the cluster

You need to install miniconda if you do not already have Anaconda Python installed in your cluster home dir. Log in to the cluster and run this command to download a miniconda installer:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  
Then this command to download and install miniconda:

    bash Miniconda3-latest-Linux-x86_64.sh

Say yes when it asks if it should run `conda init` for you.

## Creating an environment on the cluster

Log in to the cluster and run this command to create a conda envionment for your project on the cluster. This environment should contain the packages that you need for your project. Such packages may inculde:

- **grid workflow:** [gwf](https://docs.gwf.app/)
- **Jupyter:** [jupyter](https://jupyter.org/) [jupyterlab](https://jupyter.org/) jupyter_contrib_nbextensions rise nbconvert (jupyter)
- **Vectors and data frames:** numpy pandas 
- **Plotting**: matplotlib seaborn ipympl ipywidgets nodejs mpld3 plotly altair
- **Maps:** cartopy shapely fiona
- **Statistics:** scipy statsmodels 
- **Machine learning:** scikit-learn 
- **Trees:** ete3
- **Misc bioinformatics**: scikit-bio, biopython
- **Storage and indexing:** pyfaidx tabix samtools h5py
- **Graphs:** networkx
- **Gene annotation:** mygene
- **Simulation:** msprime
- **VCF files:** scikit-allel vcftools

If you run the long command below you will have access to all (and probably much more then) you need:

    conda create --name bircproject -c gwforg -c etetoolkit -c anaconda -c conda-forge -c bioconda -c kaspermunch python=3.7 slurm-jupyter jupyter jupyterlab jupyter_contrib_nbextensions rise nbconvert numpy scipy pandas h5py scikit-learn scikit-bio statsmodels matplotlib seaborn ipympl ipywidgets nodejs mpld3 plotly altair cartopy shapely fiona ete3 biopython pyfaidx networkx mygene msprime openblas scikit-allel pylint vcftools tabix samtools pip setuptools wheel twine conda-verify gwf

Should you end up needing more packages than you initially included, you can easily install them later.

**Important:** Whenever you log into the cluster to work on your project, you should activate your `bircproject` environment like this:

    conda activate bircproject
    
 When you environment is active it says `(bircproject)` on the commnad prompt instead of `(base)`.

## Visual Studio Code

If you did not do so when you installed Anaconda, you should download and instll Visual Studio Code. VS code is great for developing scripts and editing text files. Once you have installed VS code, you should install the "Remote Development" extension. You do that by clicking the funny squares in the left bar and search for "Remote Development". Once installed you can click the small green square in the lower left corner to connect to the cluster. Select "Connect current window to host" then "Add new SSH host", then type `<your_cluster_username>@login.genome.au.dk`, then select the config file `.ssh/config`. Now you can click the small green square in the lower left corner to connect to the cluster by selecting `login.genome.au.dk`. It may take a bit, but once it is done installing a remote server, you will have acces to the files in your home folder on the cluster.

## Jupyter on the cluster

Jupyter is a notebook environment where you can easily combine text, code and plots. You can read about Jupyter [here](https://jupyter.org/).

You can run a jupyter notebook in your browser from a compute node on the cluster. This way your analysis runs on the file system where your data is, and you can keep data, code and documentation in one place. [slurm-jupyter](https://github.com/kaspermunch/slurm-jupyter) is a package with a script that starts and connects to a jupyter server on compute note and forwards the web display to your local machine. 

You see the documentation on how to use and setup slurm-jupyter in the [documentation for the slurm-jupyter package](https://slurm-jupyter.readthedocs.io/en/latest/).

## Backup on the cluster

Your files on the cluster are not backed up! If you want to backup files you need to put them in a folder called BACKUP.

## Git and GitHub

Briefly what git is. Link to git and github into pages

Create a github account

Add a ssh key form the cluster

Fork this repository...

As backup

https://marklodato.github.io/visual-git-guide/index-en.html

## Backup on the cluster

link to the backup information

## Docs on cluster

https://genome.au.dk/docs/getting-started/

## Building workflows with GWF

https://gwf.app/

## Working with Jupyter notebooks

https://www.dataquest.io/blog/jupyter-notebook-tutorial/

