"""Three-column visual reports, using the supplied notebook display conventions."""
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from .core import load_fits,get_scale,CheckConfig

def plot_fits(ax,item,config):
    data,wcs=load_fits(item['file'])
    x,y=wcs.world_to_pixel_values(item['ra_deg'],item['dec_deg'])
    ny,nx=data.shape
    scale=get_scale(wcs)
    half,rp=config.plot_arcmin/scale,config.radius_arcmin/scale
    inside=np.isfinite([x,y]).all() and 0<=x<nx and 0<=y<ny
    x1,x2,y1,y2=0,nx,0,ny
    if inside:
        x1,x2=max(0,int(x-half)),min(nx,int(x+half+1))
        y1,y2=max(0,int(y-half)),min(ny,int(y+half+1))
    crop=data[y1:y2,x1:x2]
    good=crop[np.isfinite(crop)&(crop!=0)]
    lo,hi=np.percentile(good,[2,99.5]) if good.size else (0,1)
    if hi<=lo:
        hi=lo+1
    ax.imshow(crop,origin='lower',extent=(x1-.5,x2-.5,y1-.5,y2-.5),cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
    if inside:
        ax.add_patch(Circle((x,y),rp,fill=False,color='red',lw=1.5))
        ax.add_patch(Circle((x,y),8,fill=False,color='yellow',lw=1.5))
        ax.plot(x,y,'+',color='cyan',ms=9,mew=1.3)
    else:
        if np.isfinite([x,y]).all():
            sx,sy=(nx-1)/2,(ny-1)/2
            dx,dy=x-sx,y-sy
            length=np.hypot(dx,dy)
            if length:
                ax.annotate('',xy=(sx+dx/length*min(nx,ny)*.3,sy+dy/length*min(nx,ny)*.3),xytext=(sx,sy),arrowprops=dict(arrowstyle='->',color='yellow',lw=2))
        ax.text(.5,.05,'Galaxy centre outside field / invalid WCS position',color='yellow',fontsize=7,ha='center',transform=ax.transAxes)
    ax.set(xlim=(x1-.5,x2-.5),ylim=(y1-.5,y2-.5),aspect='equal',xticks=[],yticks=[])
    ax.set_title(f'{item["galaxy"]} | {item["band"]} | {item["status"]}\n{Path(item["file"]).name}\nCoverage: {item["coverage_fraction"]:.1%}',fontsize=7)

def create_pdf(table,statuses,filename,config=None):
    config=config or CheckConfig()
    subset=table.loc[table.status.isin(statuses)].reset_index(drop=True)
    with PdfPages(filename) as pdf:
        for i in range(0,max(1,len(subset)),9):
            fig=Figure(figsize=(12,12))
            FigureCanvasAgg(fig)
            axes=fig.subplots(3,3).ravel()
            batch=subset.iloc[i:i+9]
            for ax,(_,item) in zip(axes,batch.iterrows()):
                try:
                    plot_fits(ax,item,config)
                except Exception as exc:
                    ax.text(.5,.5,f'{item["galaxy"]} | {item["band"]}\n{Path(item["file"]).name}\nERROR: {exc}',ha='center',va='center',fontsize=7,wrap=True,transform=ax.transAxes)
                    ax.axis('off')
            for ax in axes[len(batch):]:
                ax.axis('off')
            if subset.empty:
                axes[0].text(.5,.5,'No observations',ha='center',transform=axes[0].transAxes)
            fig.tight_layout(pad=1.2)
            pdf.savefig(fig)
            fig.clear()
