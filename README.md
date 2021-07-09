# NLSE
Repository for project NLSE at Max Planck Institut



# UPDATE 09.07.20

**The file init_variables.py is a more organized version of the code. This will be used in a final version integrating Pandas and using better ways to share data between Callbacks on Dash. You can try it using prop_matplot.py!!**




# Python Installation:




Install a python interpreter: [https://www.python.org/downloads/](https://www.python.org/downloads/) -> Development version [3.7.9](https://www.python.org/downloads/release/python-379/ ) 

check the instalation by typing on the terminal of VS code:

**`python --version`**

## Create a Virtual environment 



 Type: 

**`python  -m venv .venv`**



then: 

**` .venv\scripts\activate`**



if you get the error "Activate.ps1 is not digitally signed. You cannot run this script on the current system." then type:

**`Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process`**



Again: 

**`.venv\scripts\activate`**

Now we should **upgrade pip:**

**`python -m pip install --upgrade pip`**

The next step is to install the packages by using command line:


**`python -m pip install -r requirements.txt`**



Please create a .gitignore file and add the following lines:

**`.venv/`**

**`.vscode/`**


##  Run the App 

You can run the App just by using: 

**`python .\index.py`**

and go to the page (default): http://127.0.0.1:8050/apps/gvd




