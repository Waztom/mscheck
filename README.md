# <a name="MScheck for locating target compound masses in mass spectra"></a>**MScheck for locating target compound masses in mass spectra**

MScheck is a python package that hunts for a target compound mass + given ion mass (eg. H+, Na+).
MScheck was created to assist with the automated mass spectrum analysis of target compounds synthesised using
a high throughput approach.

# <a name="MScheck installation"></a>**MScheck installation**

MScheck relies on rdkit for generating molecule SVG images and for calculating molecular weights<br>
Installing rdkit using conda works best followed by a pip install of MScheck<br>

1. Create a conda environment with rdkit

   > `conda create -c conda-forge -n MScheck rdkit`

2. Activate the MScheck conda environment created

   > `conda activate MScheck`

3. Pip install MScheck
   > `pip install mscheck`

# <a name="MScheck use"></a>**MScheck use**

MScheck has been tested on using Agilent LCMS files (.D) as the starting file format

1. Convert .D Agilent folder to a .d MassHunter format using Agilent's ChemStation to MassHunter Translator (B.04.00)
2.

MScheck is installed via pip, <br>

MScheck is a python package that hunts for a target compound mass + given ion mass (eg. H+, Na+),<br>

**If you wish to use XChem-CAR you are advised to use the webapp:** <br>
[URL ONCE LIVE]

**Continue with this guide if you wish to download as setup XChem-CAR
for development purpouses**<br><br>

## <a name="VisualStudioCode"></a>Visual Studio Code

these instructions are designed for Visual Studio Code which can be istalled for free from: https://code.visualstudio.com/

## <a name="GitCryptKey"></a>Git-Crypt Key

Secrets required for running CAR are encypted, to unencrpyt and run you will need the key from the XChem-CAR softwear maintainer<br><br>

# <a name="RepositoryfromGitHub"></a>Clone the "xchem-car" repository from GitHub

XChem-CAR uses GitHub for version control<br>
to get started with working on CAR clone the <em>"xchem-CAR"</em> Repository from github to your device.<br><br>

### <a name="UsefulGitHubbranches"></a>Useful GitHub branches

| Branch  | Descripton                                                 | URL                                              |
| ------- | ---------------------------------------------------------- | ------------------------------------------------ |
| Main    | Most recent, stable, release                               | https://github.com/Waztom/xchem-CAR              |
| Develop | new features will be added <br> here before being released | https://github.com/Waztom/xchem-CAR/tree/Develop |

all branches can be found: https://github.com/Waztom/xchem-CAR/branches
<br><br>

# <a name="Docker"></a>Docker

## <a name="InstallDocker"></a>Install Docker

First you'll need Docker Desktop (or the relevent Docker Engine on Linux) you can find the apprpriate download like this on https://www.docker.com/get-started

## <a name="InstallDockerCompose"></a>Install Docker Compose

once Docker is installed alos install docker compose, instructions for Mac, Windows, Linux and other options are avalible: https://docs.docker.com/compose/install/

## <a name="InstallVSCodeExtention"></a>Install VS Code Extention

Docker and Docker Compose should now be installed <br>
<em>(If on Windows/Mac, start docker desktop)</em><br>

Open Visual Studio Code<br>
to check docker is running correctly open the terminal and run:<br>

> `docker --version`<br>

you should get a response similar to:

> Docker version 18.09.2, build 6247962

In Visual Studio Code open the extiensions panel (left-hand pannle or using Ctrl+Shift+X ) and search for "<em>Remote - Containers extension</em>" and click **Install**.

Once installed a box with two arrows pointing in oposit direcitons should appere in the bottom left corner of Visual studio code
<br>
<br>

# <a name="gitcrypt"></a>git-crypt

git-crypt (https://github.com/AGWA/git-crypt) is used for encypting secrets required to run CAR
you need the appropriate <em>crypt-key</em> file from the softwear maintainer.

## <a name="InstallGitCrypt"></a>Install Git-Crypt

download the compressed git-crypt package (https://www.agwa.name/projects/git-crypt/downloads/git-crypt-0.6.0.tar.gz)

Extract the package to create the directory **git-crypt-0.6.0**

in terminal:

> `$ cd git-crypt-0.6.0`<br> >`$ make`<br> >`# make install`

## <a name="UnlockingSecrets"></a>Unlocking Secrets

once Git-Crypt is installed unlock the secrets using:

> `cd [file path to XChem-CARS on your device]`<br> >`git-crypt unlock [path to git-crypt crypt-key]`

# <a name="Startsystem"></a>Start system

### <a name="LocateReopsitory"></a>Locate Reopsitory

- in terminal change directory to your coppy of the repo :

  > `cd `<em>`[local file path to xchem-CAR repoistory]`

    </em>
    or open VS Code and go to File-> Open Folder and open the repository directory<br>

### <a name="StartRemoteContainer"></a>Start Remote Container

- start Visual Studio remote container with **Ctrl + Shift + P** and type **"Remote-containers: Open folder in container"** then click on that option. <br> ensure you have the repsoitory folder [your file path/xchem-CAR] selected and choose **"Ok"**/**"Open"**
- Your container should start to build, click on the popup notificaiton at the bottom right of visual studio to view the log/progress

### <a name="TimetoLaunch"></a>Time to Launch

- Open a new terminal that you can interact with. if the terminal is visible at the bottom of the screen click on the plus "create new intergreated terminal" or use the keybord shortcut "**Ctrl+Shift+`**" button or use the ajacent "split terminal" (or "**Ctrl+Shift+5**") button to see the new teminal ajacent to the current teminal
- you should now be in the container running Debian Linux
- in the new teminal type:

  > `cd CAR` <br> >`python3 manage.py makemigrations backend` <br> >`python3 manage.py migrate backend` <br> >`python3 manage.py runserver`<br>

  before you open the application **you must** compleate the next step, starting Celery

### <a name="StartingCelery"></a>Starting Celery

- open a new teminal the same way as last time ([see Time to Launch](#TimeToLaunch))
- in the new teminal type:
  > `celery -A CAR worker -l info`

# Opening the application

at the end of the step "[Time to Launch](#TimeToLaunch)" a address to use the visual interface should have been displayed ("http://127.0.0.1:8000/"), Ctrl+Click on the link in teminal or coppy and paste the link into your web browser to use the CAR interface

> http://127.0.0.1:8000/
