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
    filename = 'output.txt'  # The name of the file containing particle data
    xpos, ypos = read_particle_data(filename)

    # Delete every nth frame (example: every 3rd frame)
    n = 100
    xpos, ypos = keep_every_nth_frame(xpos, ypos, n)

    fig, ax = plt.subplots()
    ax.set_xlim((min(xpos) * 1.1), (max(xpos) * 1.1))
    ax.set_ylim((min(ypos) * 1.1), (max(ypos) * 1.1))
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title('Particle Trajectory')

    line, = ax.plot([], [], lw=2)

    ani = animation.FuncAnimation(fig, animate, fargs=(xpos, ypos, line),
                                  frames=len(xpos), interval=1, blit=True)

    plt.show()

if __name__ == '__main__':
    main()
