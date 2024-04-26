from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from utils import *

det_height = 600
det_radius = 200
det_center = [0,0,0]

def det_area():
    x = np.linspace(-det_radius, det_radius, 100)
    z = np.linspace(-det_height/2, det_height/2, 100)
    X, Z = np.meshgrid(x, z)
    Y = np.sqrt(det_radius**2 - X**2)
    ax.plot_surface(X+det_center[0], Y+det_center[1], Z+det_center[2], color='orange', alpha=0.3)
    ax.plot_surface(X+det_center[0], -Y+det_center[1], Z+det_center[2], color='orange', alpha=0.3)

if __name__ == '__main__':
    if not os.path.exists('./indet_body.csv'):
        poor_guys = pd.read_csv('./whole_body.csv')
        # poor_guys = poor_guys.loc[poor_guys.x**2 + poor_guys.y**2 < 100**2]
        poor_guys.set_index('muonid', inplace=True)
        poor_guys['nx'] = poor_guys['px'] / np.sqrt(poor_guys['px']**2 + poor_guys['py']**2 + poor_guys['pz']**2)
        poor_guys['ny'] = poor_guys['py'] / np.sqrt(poor_guys['px']**2 + poor_guys['py']**2 + poor_guys['pz']**2)
        poor_guys['nz'] = poor_guys['pz'] / np.sqrt(poor_guys['px']**2 + poor_guys['py']**2 + poor_guys['pz']**2)

        #find the muons that hit the detector
        for index, track in poor_guys.iterrows():
            xy_norm = np.linalg.norm(np.array([track['nx'],track['ny'],0]))
            z_test = np.linspace(500-det_height/2, 500+det_height/2, 600)
            x_test = track['x'] - track['nx']*z_test/track['nz']
            y_test = track['y'] - track['ny']*z_test/track['nz']
            for i in range(len(z_test)):
                if x_test[i]**2 + y_test[i]**2 < det_radius**2:
                    break
                if i == len(z_test)-1:
                    poor_guys.drop(index, inplace=True)
                    break

        print(len(poor_guys))
        poor_guys.to_csv('./indet_body.csv')
    else:
        poor_guys = pd.read_csv('./indet_body.csv')
        poor_guys.set_index('muonid', inplace=True)

    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(poor_guys['x'], poor_guys['y'],500, c=-poor_guys['nz'], cmap='viridis',vmin=0.0,vmax=1,s=0.5)
    ax.quiver(poor_guys['x'], poor_guys['y'],500, 200*poor_guys['nx'], 200*poor_guys['ny'], 200*poor_guys['nz'],
               colors=plt.cm.viridis(-poor_guys['nz']),arrow_length_ratio=0.1)
    det_area()
    ax.set_zlim(-300,300)
    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)
    plt.colorbar(scatter)
    
    plt.show()