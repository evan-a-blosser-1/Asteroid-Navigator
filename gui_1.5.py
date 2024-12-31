"""_summary_

Raises:
    ValueError: _description_

Returns:
    _type_: _description_
    
    
    TO-DO:

    - create hamiltonian energy plot, maybe in second tab 
    
    - make approach tab, hohmman transfer fro Earth to asteroid
        (Using jpl horizons)
        
    -
    
    
"""
import os
import sys
import numpy as np
import pandas as pd
import datetime
import trimesh
import threading
import time
from scipy.integrate import solve_ivp
from PySide6.QtCore import Qt, QDir
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                               QWidget, QLabel, QLineEdit, QPushButton, QTextEdit, QTabWidget,
                               QMenuBar, QComboBox, QScrollArea, QListView, QFileDialog, 
                               QFileSystemModel,QGridLayout,QMessageBox)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
##########################################
######################### Personal Modules
sys.dont_write_bytecode = True
import constants as C
from simulation import run_simulation
from EqPt_Stability_Sim import run_eqpt, run_eqpt_simulation
from Hohm_Horizens import run_hohmann
import Aster_Database as AD 
from MASCON_Modeler import MASCON_Modeler, plot_modeler_results

##########################################

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig, self.axis = plt.subplots(subplot_kw={"projection": "3d"})
        super().__init__(self.fig)
        self.setParent(parent)
        self.set_canvas_configuration()

    def set_canvas_configuration(self, theme='Dark Mode'):
        self.fig.clf()
        self.axis = self.fig.add_subplot(projection="3d")
        self.axis.set_xlabel('X (km)')
        self.axis.set_ylabel('Y (km)')
        self.axis.set_zlabel('Z (km)')
        self.axis.margins(0.42)
        Background = "#000000"
        Grid_Color = '#1A85FF'
        self.fig.set_facecolor(Background)
        self.axis.set_facecolor(Background)
        self.axis.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.tick_params(axis='x', colors=Grid_Color)
        self.axis.tick_params(axis='y', colors=Grid_Color)
        self.axis.tick_params(axis='z', colors=Grid_Color)
        self.axis.yaxis.label.set_color(Grid_Color)
        self.axis.xaxis.label.set_color(Grid_Color)
        self.axis.zaxis.label.set_color(Grid_Color)
        self.axis.xaxis.line.set_color(Grid_Color)
        self.axis.yaxis.line.set_color(Grid_Color)
        self.axis.zaxis.line.set_color(Grid_Color)
        plt.rcParams['grid.color'] = Grid_Color
        
        
    ###
    def plot_sphere(self, radius=1, center=(0,0,0), color='blue', alpha=0.3, label='Sun'):
        # Create sphere coordinates
        phi = np.linspace(0, 2*np.pi, 100)
        theta = np.linspace(0, np.pi, 100)
        phi, theta = np.meshgrid(phi, theta)
        
        # Convert to Cartesian coordinates
        x = radius * np.sin(theta) * np.cos(phi) + center[0]
        y = radius * np.sin(theta) * np.sin(phi) + center[1]
        z = radius * np.cos(theta) + center[2]
        
        # Plot the sphere
        self.axis.plot_surface(x, y, z, color=color, alpha=alpha, label=label)
        # Calculate orbit radius from center point
        orbit_radius = np.sqrt(center[0]**2 + center[1]**2 + center[2]**2)
        
        if orbit_radius > 0:  # Only plot ring if not at origin
            # Create ring coordinates
            t = np.linspace(0, 2*np.pi, 100)
            ring_x = orbit_radius * np.cos(t)
            ring_y = orbit_radius * np.sin(t)
            ring_z = np.zeros_like(t)
            
            # Plot the ring
            self.axis.plot(ring_x, ring_y, ring_z, color='white', 
                            alpha=0.5, linestyle='--',linewidth=0.3)
            

    def slrsys_canvas_configuration(self, theme='Dark Mode'):
        self.fig.clf()
        self.axis = self.fig.add_subplot(projection="3d")
        self.plot_sphere(radius=3, center=(0, 0, 0), color='Yellow', alpha=0.9, label='Sun')
        self.plot_sphere(radius=1, center=(5, 5, 0), color='sandybrown', alpha=0.9, label='Mercury')
        self.plot_sphere(radius=1, center=(12, 12, 0), color='rosybrown', alpha=0.9, label='Venus')
        self.plot_sphere(radius=1, center=(15, 15, 0), color='blue', alpha=0.9, label='Earth')
        self.plot_sphere(radius=1, center=(20, 20, 0), color='red', alpha=0.9, label='Mars')
        self.plot_sphere(radius=2, center=(25, 25, 0), color='orange', alpha=0.9, label='Jupiter')
        self.plot_sphere(radius=1, center=(30, 30, 0), color='sienna', alpha=0.9, label='Saturn')
        self.plot_sphere(radius=1, center=(35, 35, 0), color='cyan', alpha=0.9, label='Uranus')
        self.plot_sphere(radius=1, center=(40, 40, 0), color='purple', alpha=0.9, label='Neptune')
        self.plot_sphere(radius=1, center=(45, 45, 0), color='pink', alpha=0.9, label='Pluto')



        self.axis.set_xlabel('X (km)')
        self.axis.set_ylabel('Y (km)')
        self.axis.set_zlabel('Z (km)')
        self.axis.margins(0.42)
        Background = "#000000"
        Grid_Color = '#1A85FF'
        self.fig.set_facecolor(Background)
        self.axis.set_facecolor(Background)
        self.axis.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        self.axis.tick_params(axis='x', colors=Grid_Color)
        self.axis.tick_params(axis='y', colors=Grid_Color)
        self.axis.tick_params(axis='z', colors=Grid_Color)
        self.axis.yaxis.label.set_color(Grid_Color)
        self.axis.xaxis.label.set_color(Grid_Color)
        self.axis.zaxis.label.set_color(Grid_Color)
        self.axis.xaxis.line.set_color(Grid_Color)
        self.axis.yaxis.line.set_color(Grid_Color)
        self.axis.zaxis.line.set_color(Grid_Color)
        self.axis.set_aspect('equal', 'box')
        plt.rcParams['grid.color'] = Grid_Color


    def draw_blank(self):
        self.set_canvas_configuration()
        self.draw()

    def draw_blank_slrsystem(self):
        self.slrsys_canvas_configuration()
        self.draw()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # Title snf main window
        self.setWindowTitle("Asteroid Navigator \u2609-\u263F-\u2640-\u2295\u263D-\u2642-\u2643-\u2644-\u26E2-\u2646-\u2647")
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        ##############################
        #################### Set tabs 
        self.simulation_tab = QWidget()
        self.eqpt_tab = QWidget()
        self.modeler_tab = QWidget()
        self.database_tab = QWidget() 
        self.hohmann_tab = QWidget()
        self.mission_output_tab = QWidget()

        self.tabs.addTab(self.simulation_tab, "Simulation")
        self.tabs.addTab(self.eqpt_tab, "Equilibrium Points")
        self.tabs.addTab(self.modeler_tab, "Modeler")
        self.tabs.addTab(self.database_tab, "Database") 
        self.tabs.addTab(self.hohmann_tab, "Hohmann Transfer")
        self.tabs.addTab(self.mission_output_tab, "Mission Output")

        self.setup_simulation_tab()
        self.setup_eqpt_tab()
        self.setup_modeler_tab()
        self.setup_database_tab() 
        self.setup_hohmann_tab() 
        self.setup_mission_output_tab()
        ###############################
        ###############################
        
        ###############
        ## Main Menu ##
        ###############
        # Define menu bar
        self.menu = self.menuBar()
        # Add File 
        self.menu_file = self.menu.addMenu("File")
        exit = QAction("Exit", self, triggered=qApp.quit)
        self.menu_file.addAction(exit)
        ##################
        # About Software #
        ##################
        # Define about section of menu
        self.menu_about = self.menu.addMenu("&About")
        # Define about message box
        def About_Message():
            msg_box = QMessageBox()
            msg_box.setWindowTitle("About")
            
            html_content = """
            <div style='font-family: Arial; text-align: center; min-width: 500px;'>
            <!-- Title with Shadow -->
            <h1 style='color: #2c3e50; text-shadow: 2px 2px 4px rgba(0,0,0,0.2);'>
                Asteroid Navigator
            </h1>
            <!-- Animated Version Number -->
            <h3 style='color: #34495e; transition: all 0.3s ease;'>
                Version 1.5
            </h3>

            <!-- Styled Box with Border -->
            <div style='
                margin: 20px;
                padding: 15px;
                border: 2px solid #e1e8e3;
                border-radius: 10px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                background-color: rgba(255,255,255,0.9);
            '>
                <!-- Text with Gradient -->
                <p style='
                    background: linear-gradient(45deg, #2c3e50, #3498db);
                    -webkit-background-clip: text;
                    -webkit-text-fill-color: transparent;
                    font-weight: bold;
                '>
                    This software was made as a teaching tool for Astrodynamics!😊
                </p>

                <!-- Dedication Section -->
                <div style='
                    margin-top: 20px;
                    padding: 15px;
                    background-color: #000000;
                    border-radius: 10px;
                '>
                    <p style='color: #FFD700;'>
                        <span style='
                            font-size: 8.2em;
                            font-weight: bold;
                            border-bottom: 5px solid #FFD700;
                            padding-bottom: 15px;
                        '>Dedicated To:</span><br><br>
                        <span style='
                            font-style: italic;
                            font-size: 4.2em;
                            line-height: 2.5;
                        '>
                            my grandmother Carol Massad-Clements,<br>
                            my Uncle Patrick Clements,<br>
                            &<br>
                            Sweet Pea, your mother and I loved you very much.<br>
                            You never got to see this world, but the<br>
                            infinite cosmos is your home now.<br>
                            I hope one day I can be a father you would<br>
                            have been proud to have.<br>
                            See ya'll again someday ❤<br>
                            🥀
                        </span>
                    </p>
                </div>
            </div>

            <!-- Footer with Hover Effect -->
            <div style='
                margin-top: 20px;
                color: #95a5a6;
                transition: color 0.3s ease;
            '>
                <p style='font-size: 0.9em;'>© 2024</p>
            </div>
        </div>
            """
            
            msg_box.setText(html_content)
            msg_box.setStyleSheet("QMessageBox { min-width: 600px; }")
            msg_box.exec()
        # Assign about in menu
        about = QAction("About", self, shortcut=QKeySequence(QKeySequence.HelpContents),
                        triggered=About_Message)
        self.menu_about.addAction(about)
        
        ###############
        # MIT License #
        ###############
        # Define License message box
        def MIT_License():
            # Call message file
            with open("Other_Info/MIT_License_Evan.md", "r") as file:
                Output_Message = file.read()
                msg_box = QMessageBox()
                msg_box.setWindowTitle("MIT License")
                msg_box.setText(Output_Message)
                msg_box.exec()
        # Assign license in menu
        License = QAction("License", self, shortcut=QKeySequence(QKeySequence.HelpContents),
                        triggered=MIT_License)
        self.menu_about.addAction(License)
        
        
        
#%% Simulation Tab
###################################################
########################### Simulation Tab
    def setup_simulation_tab(self):
        self.simulation_main_layout = QHBoxLayout()

        self.simulation_control_layout = QVBoxLayout()
        


        self.asteroid_label = QLabel("Asteroid Name:")
        # self.combo_folder = QComboBox(self)
        asteroid_folder = 'Asteroids'
        def refresh_folders():
            self.combo_folder.clear()
            folders = [f for f in os.listdir(asteroid_folder) if os.path.isdir(os.path.join(asteroid_folder, f))]
            self.combo_folder.addItems(folders)
        # Override the showEvent method to refresh folders when the combo box is shown
        class RefreshableComboBox(QComboBox):
            def showEvent(self, event):
                refresh_folders()
                super().showEvent(event)
        self.combo_folder = RefreshableComboBox(self)
        self.simulation_control_layout.addWidget(self.asteroid_label)
        self.simulation_control_layout.addWidget(self.combo_folder)




        self.y0_label = QLabel("Initial y0:")
        self.y0_input = QLineEdit("0.5")
        self.simulation_control_layout.addWidget(self.y0_label)
        self.simulation_control_layout.addWidget(self.y0_input)

        self.h0_label = QLabel("Initial H0:")
        self.h0_input = QLineEdit("1.6e-9")
        self.simulation_control_layout.addWidget(self.h0_label)
        self.simulation_control_layout.addWidget(self.h0_input)

        self.days_label = QLabel("Number of Days:")
        self.days_input = QLineEdit("1.0")
        self.simulation_control_layout.addWidget(self.days_label)
        self.simulation_control_layout.addWidget(self.days_input)

        self.run_button = QPushButton("Run Simulation")
        self.run_button.clicked.connect(self.call_simulation)
        self.simulation_control_layout.addWidget(self.run_button)

        self.output = QTextEdit()
        self.simulation_control_layout.addWidget(self.output)

        self.orbit_output = QTextEdit()
        self.simulation_control_layout.addWidget(self.orbit_output)

        self.simulation_plot_canvas = PlotCanvas(self)
        self.simulation_toolbar = NavigationToolbar(self.simulation_plot_canvas, self)

        self.simulation_plot_layout = QVBoxLayout()
        self.simulation_plot_layout.addWidget(self.simulation_toolbar)
        self.simulation_plot_layout.addWidget(self.simulation_plot_canvas)

        self.simulation_main_layout.addLayout(self.simulation_control_layout)
        self.simulation_main_layout.addLayout(self.simulation_plot_layout)

        self.simulation_tab.setLayout(self.simulation_main_layout)

        self.simulation_plot_canvas.draw_blank()
        
        
    def call_simulation(self):
        folder_name = self.combo_folder.currentText()
        file_name = folder_name  # Assuming the file name is the same as the folder name
        selected_extension = "CM"  # Replace with actual extension logic

        data_file = f'Asteroids/{folder_name}/{file_name}_{selected_extension}.in'
        obj_file = f'Asteroids/{folder_name}/{file_name}.obj'
        mu_file = f'Asteroids/{folder_name}/{file_name}_mu.in'
        Const_file = f'Asteroids/{folder_name}/{file_name}_const.in'

        # Check if data_file and mu_file exist
        if not os.path.exists(data_file) or not os.path.exists(mu_file):
            self.output.append("Data file or mu file does not exist. Please create the model in the modeler tab.")
            self.run_button.setEnabled(True)
            return

        # Disable the run button to prevent multiple clicks
        self.run_button.setEnabled(False)
        self.output.append("Starting simulation...")

        # Create a thread to run the simulation
        thread = threading.Thread(target=self.run_simulation_thread, args=(data_file, obj_file, mu_file, file_name,Const_file))
        thread.start()

    def run_simulation_thread(self, data_file, obj_file, mu_file, file_name, Const_file):
        y0 = float(self.y0_input.text())
        H0 = float(self.h0_input.text())
        days = float(self.days_input.text())
        
        run_simulation(self, data_file, obj_file, mu_file, file_name, Const_file,y0,H0,days)
        # Re-enable the run button after the simulation is complete
        self.run_button.setEnabled(True)
        
#%% Equilibrium Points Tab
###################################################
########################### Equilibrium Points Tab
    def setup_eqpt_tab(self):
        self.eqpt_main_layout = QHBoxLayout()
        self.eqpt_control_layout = QVBoxLayout()



        self.eqpt_label = QLabel("Asteroid Name:")
        asteroid_folder = 'Asteroids'
        def refresh_eqpt_folders():
            self.eqpt_combo_folder.clear()
            folders = [f for f in os.listdir(asteroid_folder) if os.path.isdir(os.path.join(asteroid_folder, f))]
            self.eqpt_combo_folder.addItems(folders)      
        class RefreshableEqptComboBox(QComboBox):
            def showEvent(self, event):
                refresh_eqpt_folders()
                super().showEvent(event)
        self.eqpt_combo_folder = RefreshableEqptComboBox(self)
        self.eqpt_control_layout.addWidget(self.eqpt_label)
        self.eqpt_control_layout.addWidget(self.eqpt_combo_folder)

        ##################
        # Eq point controls 
        self.EQ_tol_label = QLabel("Tolerance:")
        self.EQ_input_tol = QLineEdit("1e-09")
        self.EQ_step_label = QLabel("Step Size:")
        self.EQ_input_step = QLineEdit("0.01")
        self.EQ_max_iter_label = QLabel("Max Iterations:")
        self.EQ_input_max_iter = QLineEdit("100000")
        self.EQ_gmesh_label = QLabel("Mesh Guess Size:")
        self.EQ_input_gmesh = QLineEdit("2")
        
        EQ_sett_grid_layout = QGridLayout()
        EQ_sett_grid_layout.addWidget(self.EQ_tol_label, 0, 0)
        EQ_sett_grid_layout.addWidget(self.EQ_input_tol, 1, 0)
        EQ_sett_grid_layout.addWidget(self.EQ_step_label, 0, 1)
        EQ_sett_grid_layout.addWidget(self.EQ_input_step, 1, 1)
        EQ_sett_grid_layout.addWidget(self.EQ_max_iter_label, 0, 2)
        EQ_sett_grid_layout.addWidget(self.EQ_input_max_iter, 1, 2)
        EQ_sett_grid_layout.addWidget(self.EQ_gmesh_label, 0, 3)
        EQ_sett_grid_layout.addWidget(self.EQ_input_gmesh, 1, 3)
        
        ########################
        # Sim Inputs 
        self.eqpt_x_label = QLabel("X (km):")
        self.eqpt_input_x = QLineEdit()
        self.eqpt_y_label = QLabel("Y (km):")
        self.eqpt_input_y = QLineEdit()
        self.eqpt_z_label = QLabel("Z (km):")
        self.eqpt_input_z = QLineEdit()
        self.eqpt_days_label = QLabel("days:")
        self.eqpt_input_days = QLineEdit()

        eqpt_grid_layout = QGridLayout()
        eqpt_grid_layout.addWidget(self.eqpt_x_label, 0, 0)
        eqpt_grid_layout.addWidget(self.eqpt_input_x, 0, 1)
        eqpt_grid_layout.addWidget(self.eqpt_y_label, 0, 2)
        eqpt_grid_layout.addWidget(self.eqpt_input_y, 0, 3)
        eqpt_grid_layout.addWidget(self.eqpt_z_label, 0, 4)
        eqpt_grid_layout.addWidget(self.eqpt_input_z, 0, 5)
        eqpt_grid_layout.addWidget(self.eqpt_days_label, 0, 6)
        eqpt_grid_layout.addWidget(self.eqpt_input_days, 0, 7)

        ################################
        ################################
        self.eqpt_control_layout.addLayout(EQ_sett_grid_layout)
        # Run button 
        self.eqpt_button = QPushButton("Run Equilibrium Point Calculation")
        self.eqpt_button.clicked.connect(self.call_eqpt_calculation)
        self.eqpt_control_layout.addWidget(self.eqpt_button)
        # Output Window
        self.eqpt_logback = QTextEdit()
        self.eqpt_logback.setReadOnly(True)
        # self.eqpt_logback.setFixedHeight(75)
        self.eqpt_control_layout.addWidget(self.eqpt_logback)

        self.eqpt_output = QTextEdit()
        # self.eqpt_output.setFixedHeight(100)
        self.eqpt_control_layout.addWidget(self.eqpt_output)

        ############################
        # Eq inputs for Linearized MASCON
        self.eqpt_control_layout.addLayout(eqpt_grid_layout)

        ############################
        # Call linearized MASCON simulation 
        self.eqpt_sim_button = QPushButton("Run Eq. point Stability Simulation")
        self.eqpt_sim_button.clicked.connect(self.call_eqpt_simulation)
        self.eqpt_control_layout.addWidget(self.eqpt_sim_button)
        ####

        ############################
        # Right Side Plot 
        self.eqpt_plot_canvas = PlotCanvas(self)
        self.eqpt_toolbar = NavigationToolbar(self.eqpt_plot_canvas, self)

        self.eqpt_plot_layout = QVBoxLayout()
        self.eqpt_plot_layout.addWidget(self.eqpt_toolbar)
        self.eqpt_plot_layout.addWidget(self.eqpt_plot_canvas)

        self.eqpt_main_layout.addLayout(self.eqpt_control_layout)
        self.eqpt_main_layout.addLayout(self.eqpt_plot_layout)

        self.eqpt_tab.setLayout(self.eqpt_main_layout)

        self.eqpt_plot_canvas.draw_blank()
        
        
        
        
    def call_eqpt_calculation(self):
        folder_name = self.eqpt_combo_folder.currentText()
        file_name = folder_name  # Assuming the file name is the same as the folder name
        selected_extension = "CM"  # Replace with actual extension logic
        data_file = f'Asteroids/{folder_name}/{file_name}_{selected_extension}.in'
        obj_file = f'Asteroids/{folder_name}/{file_name}.obj'
        mu_file = f'Asteroids/{folder_name}/{file_name}_mu.in'
        Const_file = f'Asteroids/{folder_name}/{file_name}_const.in'

        # Check if data_file and mu_file exist
        if not os.path.exists(data_file) or not os.path.exists(mu_file):
            self.eqpt_logback.append("Data file or mu file does not exist. Please create the model in the modeler tab.")
            self.eqpt_button.setEnabled(True)
            return

        # Disable the run button to prevent multiple clicks
        self.eqpt_button.setEnabled(False)
        self.eqpt_logback.append("Calculating Equilibrium Points...")

        # Create a thread to run the simulation
        thread = threading.Thread(target=self.run_eqpt_thread, args=(data_file, obj_file, mu_file, file_name,Const_file))
        thread.start()
        


    def run_eqpt_thread(self, data_file, obj_file, mu_file, file_name, Const_file):
        tol        = float(self.EQ_input_tol.text())
        max_iter   = int(self.EQ_input_max_iter.text())
        step_limit = float(self.EQ_input_step.text())
        gms        = int(self.EQ_input_gmesh.text())
        run_eqpt(self, data_file, obj_file, mu_file, file_name, Const_file, tol, max_iter, step_limit, gms)
        # Re-enable the run button after the simulation is complete
        self.eqpt_button.setEnabled(True)
        
        
        
        
    def call_eqpt_simulation(self):
        folder_name = self.eqpt_combo_folder.currentText()
        file_name = folder_name  # Assuming the file name is the same as the folder name
        selected_extension = "CM"  # Replace with actual extension logic
        data_file = f'Asteroids/{folder_name}/{file_name}_{selected_extension}.in'
        obj_file = f'Asteroids/{folder_name}/{file_name}.obj'
        mu_file = f'Asteroids/{folder_name}/{file_name}_mu.in'
        Const_file = f'Asteroids/{folder_name}/{file_name}_const.in'

        # Check if data_file and mu_file exist
        if not os.path.exists(data_file) or not os.path.exists(mu_file):
            self.eqpt_logback.append("Data file or mu file does not exist. Please create the model in the modeler tab.")
            self.eqpt_sim_button.setEnabled(True)
            return

        # Disable the run button to prevent multiple clicks
        self.eqpt_sim_button.setEnabled(False)
        self.eqpt_logback.append("simulating orbit with equilibrium points...")

        # Create a thread to run the simulation
        thread = threading.Thread(target=self.run_eqpt_sim_thread, args=(data_file, obj_file, mu_file, file_name,Const_file))
        thread.start()
        


    def run_eqpt_sim_thread(self, data_file, obj_file, mu_file, file_name, Const_file):
        x_eqpt = float(self.eqpt_input_x.text())
        y_eqpt = float(self.eqpt_input_x.text())
        z_eqpt = float(self.eqpt_input_x.text())
        run_eqpt_simulation(self, data_file, obj_file, mu_file, file_name, Const_file, x_eqpt, y_eqpt, z_eqpt)
        # Re-enable the run button after the simulation is complete
        self.eqpt_sim_button.setEnabled(True)
        
        
        
        
#%% Modeler Tab
###################################################
########################### Modeler Tab
    def setup_modeler_tab(self):
        self.modeler_main_layout = QHBoxLayout()

        self.modeler_control_layout = QVBoxLayout()



        self.modeler_asteroid_label = QLabel("Asteroid Name:")
        # self.modeler_combo_folder = QComboBox(self)
        asteroid_folder = 'Asteroids'
        def refresh_folders():
            self.modeler_combo_folder.clear()
            folders = [f for f in os.listdir(asteroid_folder) if os.path.isdir(os.path.join(asteroid_folder, f))]
            self.modeler_combo_folder.addItems(folders)
        # Override the showEvent method to refresh folders when the combo box is shown
        class RefreshableComboBox(QComboBox):
            def showEvent(self, event):
                refresh_folders()
                super().showEvent(event)
        self.modeler_combo_folder = RefreshableComboBox(self)
        self.modeler_control_layout.addWidget(self.modeler_asteroid_label)
        self.modeler_control_layout.addWidget(self.modeler_combo_folder)


        self.modeler_density_label = QLabel("Density (kg/km^3):")
        self.modeler_input_density = QLineEdit()
        self.modeler_mass_label = QLabel("Mass Accepted (kg):")
        self.modeler_input_mass = QLineEdit()
        self.modeler_volume_label = QLabel("Volume Accepted (km^3):")
        self.modeler_input_volume = QLineEdit()
        self.modeler_Spin_Rate_label = QLabel("Spin Rate (rot/hr):")
        self.modeler_input_Spin_Rate = QLineEdit("")
        self.modeler_rel_tol_label = QLabel("Relative Tolerance:")
        self.modeler_input_rel_tol = QLineEdit("1e-09")
        self.modeler_abs_tol_label = QLabel("Absolute Tolerance:")
        self.modeler_input_abs_tol = QLineEdit("1e-09")
        self.modeler_max_iter_label = QLabel("Max Iterations:")
        self.modeler_input_max_iter = QLineEdit("100000")
        self.modeler_rel_tol2_label = QLabel("Relative Tolerance 2:")
        self.modeler_input_rel_tol2 = QLineEdit("1e-14")
        self.modeler_abs_tol2_label = QLabel("Absolute Tolerance 2:")
        self.modeler_input_abs_tol2 = QLineEdit("1e-14")
        self.modeler_max_iter2_label = QLabel("Max Iterations 2:")
        self.modeler_input_max_iter2 = QLineEdit("100000")

        grid_layout = QGridLayout()
        grid_layout.addWidget(self.modeler_density_label, 0, 0)
        grid_layout.addWidget(self.modeler_input_density, 0, 1)
        grid_layout.addWidget(self.modeler_mass_label, 0, 2)
        grid_layout.addWidget(self.modeler_input_mass, 0, 3)
        grid_layout.addWidget(self.modeler_volume_label, 1, 0)
        grid_layout.addWidget(self.modeler_input_volume, 1, 1)
        grid_layout.addWidget(self.modeler_Spin_Rate_label, 1, 2)
        grid_layout.addWidget(self.modeler_input_Spin_Rate, 1, 3)
        grid_layout.addWidget(self.modeler_rel_tol_label, 2, 0)
        grid_layout.addWidget(self.modeler_input_rel_tol, 2, 1)
        grid_layout.addWidget(self.modeler_abs_tol_label, 2, 2)
        grid_layout.addWidget(self.modeler_input_abs_tol, 2, 3)
        grid_layout.addWidget(self.modeler_max_iter_label, 3, 0)
        grid_layout.addWidget(self.modeler_input_max_iter, 3, 1)
        grid_layout.addWidget(self.modeler_rel_tol2_label, 3, 2)
        grid_layout.addWidget(self.modeler_input_rel_tol2, 3, 3)
        grid_layout.addWidget(self.modeler_abs_tol2_label, 4, 0)
        grid_layout.addWidget(self.modeler_input_abs_tol2, 4, 1)
        grid_layout.addWidget(self.modeler_max_iter2_label, 4, 2)
        grid_layout.addWidget(self.modeler_input_max_iter2, 4, 3)


        self.modeler_control_layout.addLayout(grid_layout)

        self.modeler_button = QPushButton("Run Modeler Calculation")
        self.modeler_button.clicked.connect(self.call_modeler_calculation)
        self.modeler_control_layout.addWidget(self.modeler_button)

        self.running_logback = QTextEdit()
        self.running_logback.setReadOnly(True)
        self.running_logback.setFixedHeight(100)
        self.modeler_control_layout.addWidget(self.running_logback)

        self.modeler_output = QTextEdit()
        self.modeler_control_layout.addWidget(self.modeler_output)

        self.modeler_plot_canvas = PlotCanvas(self)
        self.modeler_toolbar = NavigationToolbar(self.modeler_plot_canvas, self)

        self.modeler_plot_layout = QVBoxLayout()
        self.modeler_plot_layout.addWidget(self.modeler_toolbar)
        self.modeler_plot_layout.addWidget(self.modeler_plot_canvas)

        self.modeler_main_layout.addLayout(self.modeler_control_layout)
        self.modeler_main_layout.addLayout(self.modeler_plot_layout)

        self.modeler_tab.setLayout(self.modeler_main_layout)

        self.modeler_plot_canvas.draw_blank()



    def call_modeler_calculation(self):
        folder_name = self.modeler_combo_folder.currentText()
        file_name = folder_name  # Assuming the file name is the same as the folder name

        file_2_model = f"Asteroids/{folder_name}/{file_name}.obj"

        # Check if file_2_model exists
        if not os.path.exists(file_2_model):
            self.running_logback.append("Model file does not exist. Please create the model in the modeler tab.")
            self.modeler_button.setEnabled(True)
            return

        # Disable the run button to prevent multiple clicks
        self.modeler_button.setEnabled(False)
        self.running_logback.append("Starting modeler calculation...")

        # Create a thread to run the modeler calculation
        thread = threading.Thread(target=self.run_modeler_calculation_thread, args=(file_2_model,))
        thread.start()

    def run_modeler_calculation_thread(self, file_2_model):
        self.run_modeler_calculation(file_2_model)
        # Re-enable the run button after the calculation is complete
        self.modeler_button.setEnabled(True)

    def run_modeler_calculation(self, file_2_model):
        folder_name = self.modeler_combo_folder.currentText()
        asteroid = folder_name  # Assuming the asteroid name is the same as the folder name
        Density = float(self.modeler_input_density.text())
        Mass_Accepted = float(self.modeler_input_mass.text())
        Vol_Accepted = float(self.modeler_input_volume.text())
        rel_tol = float(self.modeler_input_rel_tol.text())
        abs_tol = float(self.modeler_input_abs_tol.text())
        Max_Iter = int(self.modeler_input_max_iter.text())
        rel_tol2 = float(self.modeler_input_rel_tol2.text())
        abs_tol2 = float(self.modeler_input_abs_tol2.text())
        Max_Iter2 = int(self.modeler_input_max_iter2.text())
        Spin_Rate = float(self.modeler_input_Spin_Rate.text())

        Model_Save_Path = f"Asteroids/{folder_name}/"

        def log_callback(message):
            self.modeler_output.append(message)

        def running_log_callback(message):
            self.running_logback.append(message)

        Output_Array_M1, Gamma_Coeff, Verts, Faces = MASCON_Modeler(
            Model_Save_Path,  running_log_callback, log_callback, file_2_model,Spin_Rate, Density, Vol_Accepted, Mass_Accepted, 
            rel_tol=rel_tol, abs_tol=abs_tol, Max_Iter=Max_Iter,
            rel_tol2=rel_tol2, abs_tol2=abs_tol2, Max_Iter2=Max_Iter2)

        plot_modeler_results(self, Output_Array_M1,  Verts, Faces)

        # Re-enable the run button after the calculation is complete
        self.modeler_button.setEnabled(True)
        
#%% Database Tab
###################################################
########################### Database Tab
    def setup_database_tab(self):
        layout = QVBoxLayout()

        # Input field for file ID
        self.file_id_label = QLabel("File ID:")
        self.file_id_input = QLineEdit()
        layout.addWidget(self.file_id_label)
        layout.addWidget(self.file_id_input)

        # Input field for asteroid name
        self.asteroid_name_label = QLabel("Asteroid Name:")
        self.asteroid_name_input = QLineEdit()
        layout.addWidget(self.asteroid_name_label)
        layout.addWidget(self.asteroid_name_input)

        # Button to fetch the asteroid model
        self.fetch_button = QPushButton("Fetch Asteroid Model")
        layout.addWidget(self.fetch_button)

        # Connect the button to the fetch function
        self.fetch_button.clicked.connect(self.fetch_asteroid_model)

        # New text box for download status
        self.download_status = QTextEdit()
        self.download_status.setReadOnly(True)
        self.download_status.setFixedHeight(100)  # Set a fixed height
        self.download_status.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)  # Ensure vertical scrollbar is always on
        layout.addWidget(self.download_status)

        # File explorer view for Asteroids folder
        self.file_model = QFileSystemModel()
        self.file_model.setRootPath(QDir.currentPath())
        self.file_view = QListView()
        self.file_view.setModel(self.file_model)
        self.file_view.setRootIndex(self.file_model.index("Asteroids"))
        layout.addWidget(self.file_view)

        # Button to run the database query
        self.database_button = QPushButton("Run Database Query")
        layout.addWidget(self.database_button)

        # Scroll area for query results
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        layout.addWidget(self.scroll_area)

        # Output field for query results
        self.database_output = QTextEdit()
        self.database_output.setReadOnly(True)
        self.scroll_area.setWidget(self.database_output)

        # Set the layout for the tab
        self.database_tab.setLayout(layout)

        # Connect the button to the query function
        self.database_button.clicked.connect(self.run_database_query)

    def fetch_asteroid_model(self):
        asteroid_id = self.file_id_input.text()
        asteroid_name = self.asteroid_name_input.text()
        if asteroid_id and asteroid_name:
            # Disable the fetch button to prevent multiple clicks
            self.fetch_button.setEnabled(False)
            self.download_status.append("Starting download, please wait for files...")

            # Create a thread to run the download
            thread = threading.Thread(target=self.run_fetch_asteroid_model_thread, args=(asteroid_id, asteroid_name))
            thread.start()
        else:
            self.database_output.setPlainText("Please enter both File ID and Asteroid Name.")

    def run_fetch_asteroid_model_thread(self, asteroid_id, asteroid_name):
        def log_callback(message):
            self.download_status.append(message)

        try:
            AD.download_obj_file(asteroid_id, asteroid_name, log_callback)
            log_callback("Download completed successfully.")
        except Exception as e:
            log_callback(f"Error: {str(e)}")
        finally:
            # Re-enable the fetch button after the download is complete
            self.fetch_button.setEnabled(True)


    def run_database_query(self):
        df_asteroids = AD.list_asteroids()
        
        # Organize the DataFrame by category
        categories = df_asteroids['Category'].unique()
        output_text = ""
        
        for category in sorted(categories):
            df_category = df_asteroids[df_asteroids['Category'] == category]
            output_text += f"Category {category}:\n"
            output_text += df_category[['ID', 'Name']].to_string(index=False)
            output_text += "\n\n"
        
        self.database_output.setPlainText(output_text)
        
        
        
#%% Hohmann Transfer Tab

###################################################
########################### Hohmann Transfer Tab
    def setup_hohmann_tab(self):
        self.hohmann_main_layout = QHBoxLayout()

        self.hohmann_control_layout = QVBoxLayout()
        


        self.asteroid_label = QLabel("Asteroid Name:")
        # self.hohmann_combo_folder = QComboBox(self)
        asteroid_folder = 'Asteroids'
        def refresh_folders():
            self.hohmann_combo_folder.clear()
            folders = [f for f in os.listdir(asteroid_folder) if os.path.isdir(os.path.join(asteroid_folder, f))]
            self.hohmann_combo_folder.addItems(folders)
        # Override the showEvent method to refresh folders when the combo box is shown
        class RefreshableComboBox(QComboBox):
            def showEvent(self, event):
                refresh_folders()
                super().showEvent(event)
        self.hohmann_combo_folder = RefreshableComboBox(self)
        self.hohmann_control_layout.addWidget(self.asteroid_label)
        self.hohmann_control_layout.addWidget(self.hohmann_combo_folder)




        self.hohm_y0_label = QLabel("Initial y0:")
        self.hohm_y0_input = QLineEdit("0.5")
        self.hohmann_control_layout.addWidget(self.hohm_y0_label)
        self.hohmann_control_layout.addWidget(self.hohm_y0_input)

        self.hohm_h0_label = QLabel("Initial H0:")
        self.hohm_h0_input = QLineEdit("1.6e-9")
        self.hohmann_control_layout.addWidget(self.hohm_h0_label)
        self.hohmann_control_layout.addWidget(self.hohm_h0_input)

        self.hohm_days_label = QLabel("Number of Days:")
        self.hohm_days_input = QLineEdit("1.0")
        self.hohmann_control_layout.addWidget(self.hohm_days_label)
        self.hohmann_control_layout.addWidget(self.hohm_days_input)

        self.hohmann_button  = QPushButton("Run Hohmann Transfer Calculation")
        self.hohmann_button .clicked.connect(self.call_hohmann)
        self.hohmann_control_layout.addWidget(self.hohmann_button )

        self.hohmann_output = QTextEdit()
        self.hohmann_control_layout.addWidget(self.hohmann_output)

        self.hohmann_orbit_output = QTextEdit()
        self.hohmann_control_layout.addWidget(self.hohmann_orbit_output)

        self.hohmann_plot_canvas = PlotCanvas(self)
        self.hohmann_toolbar = NavigationToolbar(self.hohmann_plot_canvas, self)

        self.hohmann_plot_layout = QVBoxLayout()
        self.hohmann_plot_layout.addWidget(self.hohmann_toolbar)
        self.hohmann_plot_layout.addWidget(self.hohmann_plot_canvas)

        self.hohmann_main_layout.addLayout(self.hohmann_control_layout)
        self.hohmann_main_layout.addLayout(self.hohmann_plot_layout)

        self.hohmann_tab.setLayout(self.hohmann_main_layout)

        self.hohmann_plot_canvas.draw_blank_slrsystem()
        
        
    def call_hohmann(self):
        folder_name = self.hohmann_combo_folder.currentText()
        file_name = folder_name  # Assuming the file name is the same as the folder name
        selected_extension = "CM"  # Replace with actual extension logic

        data_file = f'Asteroids/{folder_name}/{file_name}_{selected_extension}.in'
        obj_file = f'Asteroids/{folder_name}/{file_name}.obj'
        mu_file = f'Asteroids/{folder_name}/{file_name}_mu.in'
        Const_file = f'Asteroids/{folder_name}/{file_name}_const.in'

        # Check if data_file and mu_file exist
        if not os.path.exists(data_file) or not os.path.exists(mu_file):
            self.output.append("Data file or mu file does not exist. Please create the model in the modeler tab.")
            self.hohmann_button .setEnabled(True)
            return

        # Disable the run button to prevent multiple clicks
        self.hohmann_button .setEnabled(False)
        self.output.append("Starting simulation...")

        # Create a thread to run the simulation
        thread = threading.Thread(target=self.run_hohmann_thread, args=(data_file, obj_file, mu_file, file_name,Const_file))
        thread.start()

    def run_hohmann_thread(self, data_file, obj_file, mu_file, file_name, Const_file):
        run_hohmann(self, data_file, obj_file, mu_file, file_name, Const_file)
        # Re-enable the run button after the simulation is complete
        self.hohmann_button .setEnabled(True)
                
        
        
        
#%% Mission Output Tab
###################################################
########################### Mission output tab
    def setup_mission_output_tab(self):
        layout = QVBoxLayout()

        self.mission_output = QTextEdit()
        layout.addWidget(self.mission_output)

        self.save_button = QPushButton("Save Mission Output")
        self.save_button.clicked.connect(self.save_mission_output)
        layout.addWidget(self.save_button)

        self.preview_button = QPushButton("Preview Mission Output")
        self.preview_button.clicked.connect(self.preview_mission_output)
        layout.addWidget(self.preview_button)

        self.mission_output_tab.setLayout(layout)

    def preview_mission_output(self):
        simulation_output = self.orbit_output.toPlainText()
        modeler_output = self.modeler_output.toPlainText()
        EqPt_output = self.eqpt_output.toPlainText()
        
        # hohmann_output = self.hohmann_output.toPlainText()
        #
        # "Hohmann Transfer Output:\n" + hohmann_output + "\n\n"
        preview_text = (
            "Modeler Output:\n" + modeler_output + "\n\n" +
            "Simulation Output:\n" + simulation_output + "\n\n" +
            "Equilibrium Points Output:\n" + EqPt_output + "\n\n"
        )

        self.mission_output.setPlainText(preview_text)

    def save_mission_output(self):
        mission_output_text = self.mission_output.toPlainText()

        with open("mission_output.txt", "w") as file:
            file.write(mission_output_text)


#%% Main Function
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())