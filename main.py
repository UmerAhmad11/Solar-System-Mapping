import tkinter as tk
from tkinter import ttk, messagebox
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from skyfield.api import load
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # Import Axes3D for 3D plotting
from vpython import sphere, vector, color, scene, rate, label
import multiprocessing  # Import multiprocessing
from skyfield.data import mpc
import pandas as pd

# Load the ephemeris data using HORIZONS
eph = load('de441.bsp')  # Use a high-precision ephemeris from Skyfield
ts = load.timescale()

# Define celestial bodies
sun = eph['Sun']
earth = eph['Earth']
planets = {
    'Mercury': eph['Mercury Barycenter'],
    'Venus': eph['Venus Barycenter'],
    'Earth': eph['Earth'],
    'Moon': eph["Moon"],
    'Mars': eph['Mars Barycenter'],
    'Jupiter': eph['Jupiter Barycenter'],
    'Saturn': eph['Saturn Barycenter'],
    'Uranus': eph['Uranus Barycenter'],
    'Neptune': eph['Neptune Barycenter'],
    'Pluto': eph['Pluto Barycenter'],
}

# Load comet data
data = load.open(mpc.COMET_URL)
comets = mpc.load_comets_dataframe(data)
# Convert 'magnitude_k' to numeric values, forcing errors to NaN
comets['magnitude_k'] = pd.to_numeric(comets['magnitude_k'], errors='coerce')

comets = comets.sort_values('magnitude_k', ascending=False)


# Extract comet names and positions
bright_comets = comets.head(10) 



# Function to get live positions relative to the Sun
def get_live_positions():
    now = datetime.utcnow()  # Current UTC time
    t = ts.utc(now.year, now.month, now.day, now.hour, now.minute, now.second)

    positions = []
    coordinates = []
    
    for name, planet in planets.items():
        # Get the position in RA/Dec and distance in AU
        astrometric = sun.at(t).observe(planet)
        ra, dec, distance = astrometric.radec()  # RA, Dec, and distance
        
        # Convert to 3D Cartesian coordinates
        x = distance.au * np.cos(dec.radians) * np.cos(ra.radians)
        y = distance.au * np.cos(dec.radians) * np.sin(ra.radians)
        z = distance.au * np.sin(dec.radians)
        
        # Append the planet's 3D position and name
        coordinates.append((name, x, y, z))
        positions.append((name, distance))  # Keep distance for displaying
        
    return positions, coordinates


# Function to get a planet's position relative to Earth
def ref_earth(planet_name):
    if planet_name not in planets:
        print(f"Invalid planet name: {planet_name}")
        messagebox.showerror("Error", f"Invalid planet name: {planet_name}")
        return
    
    planet = planets[planet_name]
    now = datetime.utcnow()
    t = ts.utc(now.year, now.month, now.day, now.hour, now.minute, now.second)
    position = earth.at(t).observe(planet).apparent()
    x, y, z = position.position.au
    print(f"Position of {planet_name} relative to Earth (in AU):")
    print(f"X: {x:.2f}, Y: {y:.2f}, Z: {z:.2f}")

    return f"X: {x:.2f}, Y: {y:.2f}, Z: {z:.2f}"


# Function to update positions in the GUI and plot the solar system
def update_positions():
    positions, coordinates = get_live_positions()
    
    # Update text box with live positions
    text_box.delete(1.0, tk.END)
    text_box.insert(tk.END, f"Live Positions of Celestial Bodies at {datetime.utcnow()} (in AU):\n")
    text_box.insert(tk.END, f"{'Body':<10} {'Distance (AU)':>20}\n")
    text_box.insert(tk.END, "-" * 70 + "\n")
    
    for name, distance in positions:
        text_box.insert(tk.END, f"{name:<10} {distance.au:>20.8f}\n")
    
    text_box.insert(tk.END, "\nUpdating in 1 second...\n")
    
    # Schedule the next update in 1 second
    root.after(1000, update_positions)


# Function to plot solar system in 3D using Matplotlib
def plot_solar_system(coordinates):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the Sun at the origin
    ax.scatter(0, 0, 0, color='yellow', s=200, label='Sun')  # Large yellow dot for the Sun
    
    # Plot each planet's position
    for planet, x, y, z in coordinates:
        ax.scatter(x, y, z, label=planet)  # Plot each planet as a point
        ax.text(x, y, z, planet)  # Label each planet
    
    # Labels for axes
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
    
    # Title and legend
    ax.set_title('Solar System - 3D View')
    ax.legend()

    plt.show()


# Function to display the top 10 brightest comets
def display_comets():
    comet_text_box.delete(1.0, tk.END)
    comet_text_box.insert(tk.END, "Top 10 Brightest Comets:\n")
    
    # Dynamically determine the max width for the 'Name' column
    max_name_length = max(bright_comets['designation'].apply(len))
    
    # Format header with dynamic width for 'Name'
    comet_text_box.insert(tk.END, f"{'Name':^{max_name_length}}  {'Magnitude':^10}  {'Perihelion Distance (AU)':^25}  {'Perihelion Date':^20}\n")
    comet_text_box.insert(tk.END, "-" * (max_name_length + 60) + "\n")
    
    # Loop through the comets to insert the data
    for _, row in bright_comets.iterrows():
        name = row['designation']
        magnitude = row['magnitude_k']
        perihelion_distance = row['perihelion_distance_au']  # Extract perihelion distance
        perihelion_year = int(row['perihelion_year'])
        perihelion_month = int(row['perihelion_month'])
        perihelion_day = int(row['perihelion_day'])

        # Format the date properly
        perihelion_date = datetime(perihelion_year, perihelion_month, perihelion_day).strftime('%Y-%m-%d')

        # Format the data into fixed-width, centered columns
        comet_text_box.insert(tk.END, f"{name:^{max_name_length}}  {magnitude:^10.2f}  {perihelion_distance:^25.2f}  {perihelion_date:^20}\n")

# Function to create 3D solar system simulation using VPython
# Function to create 3D solar system simulation using VPython
def create_vpython_solar_system():
    scene.width = 1480
    scene.height = 780
    scene.title = "3D Solar System Simulation"
    scene.camera.pos = vector(1, 1, 1)  # Starting position of the camera
    scene.camera.axis = vector(-1, -1, -1)  # Direction the camera is pointing

    # Create the Sun at the origin
    sun_sphere = sphere(pos=vector(0, 0, 0), radius=2.18, color=color.yellow, emissive=True, make_trail=False)
    
    # Define colors for the planets
    planet_colors = {
        'Mercury': color.gray(0.5),
        'Venus': color.yellow,
        'Earth': color.blue,
        'Mars': color.red,
        'Jupiter': color.orange,
        'Saturn': color.yellow,
        'Uranus': color.cyan,
        'Neptune': color.blue,
        'Pluto': color.white,
        'Moon': color.white
    }
    
    # Create spheres for each planet
    planet_spheres = {}
    planet_labels = {}

    # Scaling factor to adjust distances (scaled by a factor of 10 for better visualization)
    scale_factor = 10  # You can adjust this factor to scale the positions

    # Loop over planets to create a sphere for each one at the correct position
    positions, coordinates = get_live_positions()
    for name, x, y, z in coordinates:
        # Scale the position
        x *= scale_factor
        y *= scale_factor
        z *= scale_factor
        
        if name == 'Moon':
            # Make the Moon smaller than the planets
            radius = 0.005  # Moon's radius is smaller than planets
        elif name == 'Earth':
            radius = 0.02  # Default radius for planets
        elif name == 'Mercury':
            radius = 0.007
        elif name == 'Mars':
            radius = 0.01
        elif name == 'Jupiter':
            radius = 0.22
        elif name == 'Saturn':
            radius = 0.18
        elif name == 'Venus':
            radius = 0.02
        elif name == 'Neptune':
            radius = 0.08
        elif name == 'Uranus':
            radius = 0.08
        elif name == 'Pluto':
            radius = 0.004
        
        planet_spheres[name] = sphere(
            pos=vector(x, y, z),   # Use the scaled 3D coordinates for positioning
            radius=radius,         # Size of the planet or moon
            color=planet_colors.get(name, color.white),  # Set color for each planet
            make_trail=False        # Make a trail to show the movement
        )
    
         # Create a label for each planet
        planet_labels[name] = label(
            pos=planet_spheres[name].pos + vector(0, 0, 0.05),  # Position the label slightly above the planet
            text=name,                                            # Label text is the planet's name
            xoffset=10,                                           # Horizontal offset to make the text more readable
            yoffset=10,                                           # Vertical offset
            space=20,                                             # Space between the label and the sphere
            height=8,                                            # Font size
            border=6,                                             # Border width
            font='sans',                                          # Font type
            color=color.white                                      # Label color
        )

    while True:
        positions, coordinates = get_live_positions()
        for name, x, y, z in coordinates:
            # Scale the position
            x *= scale_factor
            y *= scale_factor
            z *= scale_factor

            planet_spheres[name].pos = vector(x, y, z)
            planet_labels[name].pos = vector(x, y, z) + vector(0, 0, 0.05)
            
            if name == 'Moon':
                earth_pos = planet_spheres['Earth'].pos
                distance = 0.09  # You can adjust this distance
                moon_offset = vector(0, distance, 0)  # Offset along the z-axis (you can modify this direction)
                planet_spheres[name].pos = earth_pos + moon_offset
                planet_labels[name].pos = planet_spheres[name].pos + vector(0, 0, 0.05)  # Update label position
            
            # Slow down the animation to a reasonable frame rate
            rate(30)  # 30 frames per second



# Function to start VPython simulation in a separate process
def start_vpython_simulation():
    simulation_process = multiprocessing.Process(target=create_vpython_solar_system)
    simulation_process.start()


# Set up Tkinter GUI
root = tk.Tk()
root.title("Solar System Visualization")
root.geometry('1000x700')

# Text box for live positions
text_box = tk.Text(root, height=15, width=75)
text_box.pack()

# Comet text box
comet_text_box = tk.Text(root, height=15, width=90)
comet_text_box.pack()

# Button to display comets
comet_button = tk.Button(root, text="Show Brightest Comets", command=display_comets)
comet_button.pack()

# Planet selection dropdown
planet_label = tk.Label(root, text="Select a Planet:")
planet_label.pack()

planet_dropdown = ttk.Combobox(root, values=list(planets.keys()))
planet_dropdown.pack()

# Button to get selected planet's position relative to Earth
planet_position_button = tk.Button(root, text="Get Position", command=lambda: get_selected_planet_position())
planet_position_button.pack()

# Label to display the selected planet's position relative to Earth
planet_position_label = tk.Label(root, text="Position will appear here.")
planet_position_label.pack()

# Function to update the selected planet's position
def get_selected_planet_position():
    selected_planet = planet_dropdown.get()
    planet_position = ref_earth(selected_planet)
    planet_position_label.config(text=f"Position of {selected_planet}: {planet_position}")

# Button to open 3D Matplotlib Solar System
plot_button = tk.Button(root, text="Plot Solar System in 3D", command=lambda: plot_solar_system(get_live_positions()[1]))
plot_button.pack()

# Button to start VPython simulation
simulation_button = tk.Button(root, text="Start VPython Simulation", command=start_vpython_simulation)
simulation_button.pack()

# Start live position updates
update_positions()

# Start the GUI event loop
root.mainloop()
