from matplotlib import legend
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

from LocalEnvAL import const



def single_atom_rdf_plot(r, data: np.ndarray, title: str = 'unknown'):
    
    fig = go.Figure()
    r = r[:np.shape(data)[0]]
    for i, key in zip(range(data.shape[1]), const.PerTab.keys()):
        x = r
        y = np.zeros(data.shape[0]) + i 
        z = data[:,i]
        is_show_legend = True if np.max(z > 1e-6) else False
        alpha = 0.8  if is_show_legend else 0.1
        

        fig.add_trace(
            go.Scatter3d(
            x=x, y=y, z=z,
            mode= 'lines',
            line= dict(width=3),
            name= key,
            opacity= alpha,
            showlegend=is_show_legend)
        )
    
    fig.update_layout(
        title={
            'text': title,
            'y': 0.85,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        legend =dict(
            x=0.1,  # x坐标（0到1之间）
            y=1,    # y坐标（0到1之间）
            xanchor='left',  # x锚点
            yanchor='top'      # y锚点
        ),
        scene=dict(
            xaxis_title='r(Å)',
            yaxis_title='Z number',
            zaxis_title=r'$\rho$'
        ),
        # width=720,
        # height=480
    )
    
    return fig

def rcut_plot(r, rcut, rdf_fig:go.Figure, z_max:float=5):
    fig = rdf_fig
    y = np.linspace(0, len(const.PerTab.keys()), 10)
    z = np.linspace(0, z_max, 10)
    Y, Z = np.meshgrid(y, z)
    x = np.full_like(Y, rcut)
    fig.add_trace(go.Surface(
        x=x, y=Y, z=Z,
        opacity=0.5,
        colorscale= 'Blues',
        showscale=False,
        name="r_cut",
        showlegend=True
    ))
    
    
    return fig
    