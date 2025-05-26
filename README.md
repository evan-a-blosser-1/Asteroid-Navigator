
![image](https://github.com/user-attachments/assets/ed311138-f27b-4c6b-83ad-2ba039f95040)

---
In Early Release!
---

This is a continued project from undergraduate to graduate research in astrodynamics. It serves as a mission design aid and a teaching guide for those eager to learn! :)

# Mass Concentration (MASCON) Models for Asteroids

The method of using models to approximate gravitational fields of distant small objects like asteroids is essential for space mission design. Usually the point approximation can be used at a distant, yet once we journey to these asteroids we must account for perturbations (or acting forces) caused by the asteroids shape and size!



![image](https://github.com/user-attachments/assets/a99f5fb4-0f7a-4638-a44d-8fe61ee57029)



# Hamiltonian Mechanics:

The particles are set in motion utilizing Hamiltonian Mechanics:

![image](https://github.com/user-attachments/assets/ed23e36b-9e2f-4e46-a19b-0d6a743eda08)

Near small-bodies like asteroids, classical orbital elements fail, making the ideal choice for the mechanics to be Hamiltonian. Using the following equation for Hamiltonian energy we solve for the velocity of our particle using the assumption that we are setting off from the y-axis in the positive x-direction with a positive $\dot{x}$ velocity. Thus we assume that  $\dot{y}$ = $\dot{z}$ = x = z = 0:


![image](https://github.com/user-attachments/assets/abdec320-d568-46ef-82a8-aca8a31d56ff)

From this the value for $\dot{x}$ is found and used as the initial velocity of the particle. 

# How To Use Software:

The sections of this software are the Database, Modeler, Simulation, Equilibrium Points, and Mission Output. (Hohmann Transfer Not Complete)

## Database 

By clicking the Run Database Query, the software will fetch asteroid shape models from a known online database, given internet access. From this the ID and Asteroid name should be typed into the input fields by the user before using the Fetch Asteroid Model button to begin downloading all shape models available. 

![image](https://github.com/user-attachments/assets/6f60b9be-b293-4335-b3ec-dffbc6794fa0)


## Modeler 

This requires a little research from the user, in which the various parameters for the asteroid must be found and used to create a MASCON model of the asteroid. These values may not be widely known, found in research papers or online, and may only be best approximations. 

Finding the *accepted density & mass* of the asteroid along, with the *accepted volume*, the model can be created. There is also a place to enter the rotation rate of the asteroid which will be used in simulations. 

Other values control the precision of calculations and are preset. However, the user can adjust to correct errors in the model that can arise if the model is a large file or very irregular in shape.  

![image](https://github.com/user-attachments/assets/eb218506-e033-406b-b664-85e21340ef44)


## Simulating Orbits

Simulations utilizing Hamiltonian mechanics can be conducted once the model is created. Simply selected a desired initial location on the y-axis in kilometers, and hamiltonian energy in km/s:

![image](https://github.com/user-attachments/assets/fd6ecf3f-ddef-444b-8af3-439119e9e017)



## Equilibrium Points

This section allows for the finding of equilibrium points around the asteroid using the Newton-Raphson method in 3-Dimensions. The controls allow for refinement of the solution.

- Tolerance: controls the definition of the equilibrium to zero. The smaller the value the more precise the solution is.

- Step Size: controls how quickly the iterations move, set to 1 to allow to null the step size and allow for the Newton-Raphson to go at any pace. This may decrease precision. 

- Max iterations: controls how many times this will run before ending without finding a solution. 

- Mesh Guess Size: controls the mesh grid size for initial guesses. Increasing this exponentially increases the amount of solutions, but does not guarantee they will be unique! Keeping this low will increase speed and reduced redundant solutions. 



![3](https://github.com/user-attachments/assets/1f68deaf-d370-4638-855f-47361ee5f384)



## Mission Output

This section records the results of each section from the modeler's data to the simulation.

![4](https://github.com/user-attachments/assets/666255ea-e6db-4f5a-a2eb-cb49b9ccd960)


## Hohmann Transfer: Coming Soon!

![hohm](https://github.com/user-attachments/assets/54953297-c4f9-4098-b70e-057fb8121eb9)

 
---

## Installing `Asteroid Navigator`
  The `Asteroids_Navigator_setup.exe` can be downloaded from the "Releases" tab located on the right of this page.
  The current release is **Version_1.5**

  1. Download the `Asteroids_Navigator_setup.exe` file
  2. Locate the file and run as administrator
      - Windows will warn you that I am an Unknown publisher:
        
        ![image](https://github.com/evan-a-blosser-1/Asteroid_Shape_Model_GUI_1.0/assets/85218360/5bf24413-4d60-49d8-91e0-4affb43f8df8)
        
      - Simply click **more info** and then **Run anyway**:
        
        ![image](https://github.com/evan-a-blosser-1/Asteroid_Shape_Model_GUI_1.0/assets/85218360/0aae3d60-67a4-434b-8a8f-d52984f1683f)
        
3. Next just follow the setup you can choose where to install, and if you want a desktop shortcut.

   - I recommend creating the shortcut for ease of locating the software. 


## How To Run: *if you choose not to install* 
  1) Using an Interactive Development Environment (IDE)
     - download all the files, this can be done using git clone https://github.com/evan-a-blosser-1/Asteroid-Navigator.git
     - Then simply run `gui_1.5.py` 
