# OperationalResearchCW
A repo to keep together all of our coursework for module 21MAC175
## Authors
- James Taylor
- Sara Rubio Hermosa
- Oli Suttcliffe
- Lizzie Adam

## Setup Guide
This repo consists of 3 scripts.

`SimplexMatrixMethod.m`
`questions.m`
`graphs.m`

The simpexMatrixMethod is exactly what it says on the tin. its the simplex matrix method we have implemented in matlab. This script is NOT to be run. This is imported into the other two scripts and used to run them. This is for readability purposes. You can read qustion one in english, hit run and up will pop your result. You will get an error if you just run simplexMatrixMethod by itself. Think of this file like building blocks which we will use to build and produce our questions and graphs.

The questions file is where all our questions are stored. The first one which is uncommented is the question we were asked to solve. The preceeding questions are ones whoose solutions are infeasilbe/unbounded etc. Our code handles all of these. Feel free to uncomment out the questions you want to run (you can several questions at once).

The graphs file is the file which will produce graphs. Nothing needs to be changed or edited, just hit the big green run button and up should pop four graphs. You can see screen shots of these in out report.

What if something doesn't work?

- Make sure your PATH is correct. Matlab will do this automatically and ask you if you would like to add xyz.m to your PATH. Click yes to this alert. This is matlab trying to find the file to run.
- Make sure your workspace is clear. Matlab will save all variables to your workspace (this can be seen on the right hand side of the matlab editor). This is so you can run scripts and then paly around with the variables in the console. However if there are variables from previous scripts you have run which our named the same as ours (ie two groups might name the A matrix A or the b matrix b since it is the most logical thing to do) this might interfer. Our scripts SHOULD clear your work space before they are ran but if they do not simple type `clear` into your matlab console. 

If you do have any problems running this script please email us. This repo has been tested on several computers in hasselgrave all of which seem to work. 


