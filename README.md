
# DNRS Assignment-4

Project details:
https://github.com/syedjameel/Dynamic-Models-Stanford-Manipulator/blob/de44423261d49169cc829d1952a2c96be5d94286/Home%20task%204_%20Dynamics.pdf

## Main Task:
1. Euler-Lagrange Dynamic Model
2. Plot the Torques for Euler-Lagrange
3. Newton-Euler Dynamic Model
4. Plot the Torques for Newton-Euler

## Bonus Task
1. Drive the Manipulator using Trapezoidal Profile
2. Drive the Manipulator using Polynomial Profile
3. Plot the Torque graphs for both


## Description:
This Assignment consists of multiple files. The helper files are ```SimpleTranformations.py```, ```StanfordManipulatorKinematics.py```, 
```newton_euler.py```, ```trajectory_planning.py``` and ```utils.py```.
And this Assignment is run by a single driver code ```main.ipynb``` in the main folder.

Also Note that the Report is included in the ```main.ipynb``` file.


### Local Setup

* Clone project using command:
```angular2html
https://github.com/syedjameel/Dynamic-Models-Stanford-Manipulator.git
```

### Directory Structure:

```
Dynamic-models-for-Stanford-Manipulator/
├── Home task 4_ Dynamics.pdf
├── main.ipynb
├── newton_euler.py
├── README.md
├── screenshots
│   ├── euler-lagrange.png
│   ├── newton-euler1.png
│   └── newton-euler2.png
├── SimpleTranformations.py
├── StanfordManipulatorKinematics.py
├── trajectory_planning.py
└── utils.py
```


### How to Run

1. Run the following command to navigate to the directory using:

   ```shell
   cd Dynamic-Models-Stanford-Manipulator/ 
   ```
2. Run the following command to install dependent libraries:

   ```shell
   pip3 install -r requirements.txt
   ```

3. Finally, Run the ```main.ipynb``` file [In ```Jupyter notebook``` or ```colab```]
