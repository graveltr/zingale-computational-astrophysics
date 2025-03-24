import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Function to read data from the text file
def read_particle_data(filename):
    data = np.loadtxt(filename)
    xpos = data[:, 0]
    ypos = data[:, 1]
    return xpos, ypos

# Animation function
def animate(i, xpos, ypos, line):
    line.set_data(xpos[:i], ypos[:i])
    return line,

# Function to delete every nth frame
def keep_every_nth_frame(xpos, ypos, n):
    # Use list comprehension to skip every nth element
    xpos = [xpos[i] for i in range(len(xpos)) if (i + 1) % n == 0]
    ypos = [ypos[i] for i in range(len(ypos)) if (i + 1) % n == 0]
    return xpos, ypos

# Main function
def main():
    # Filenames of the particle data files
    filename1 = 'rk-output.txt'
    filename2 = 'euler-output.txt'
    
    # Read data for the two particles
    xpos1, ypos1 = read_particle_data(filename1)
    xpos2, ypos2 = read_particle_data(filename2)
    
    # Keep every nth frame (example: every 3rd frame)
    n = 100
    xpos1, ypos1 = keep_every_nth_frame(xpos1, ypos1, n)
    xpos2, ypos2 = keep_every_nth_frame(xpos2, ypos2, n)
    
    fig, (ax1, ax2) = plt.subplots(1, 2)  # Create side-by-side subplots
    
    # Set up the first subplot for particle 1
    ax1.set_xlim((min(xpos1) * 1.1), (max(xpos1) * 1.1))
    ax1.set_ylim((min(ypos1) * 1.1), (max(ypos1) * 1.1))
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position')
    ax1.set_title('Particle 1 Trajectory')
    line1, = ax1.plot([], [], lw=2)
    
    # Set up the second subplot for particle 2
    ax2.set_xlim((min(xpos2) * 1.1), (max(xpos2) * 1.1))
    ax2.set_ylim((min(ypos2) * 1.1), (max(ypos2) * 1.1))
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Y Position')
    ax2.set_title('Particle 2 Trajectory')
    line2, = ax2.plot([], [], lw=2)
    
    # Create the animation for both particles
    ani1 = animation.FuncAnimation(fig, animate, fargs=(xpos1, ypos1, line1), frames=len(xpos1), interval=10, blit=True)
    ani2 = animation.FuncAnimation(fig, animate, fargs=(xpos2, ypos2, line2), frames=len(xpos2), interval=10, blit=True)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
