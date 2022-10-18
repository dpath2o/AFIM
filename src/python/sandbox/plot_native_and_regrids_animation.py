

for i in range(90*24):
    F_png = f"./GRAPHICAL/animations/native_sequence/Python_Animation_03_frame_{i:04}.png"
    if not os.path.exists(F_png):
        print('plotting: '+str(t2m_na.coords['time'].values[i])[:13])
        t2m.isel(time=i).t2m.plot( figsize = (12,6),
                                   vmin=190, vmax=340,
                                   cmap='coolwarm',    # Change the colormap back to 'bwr'
                                   cbar_kwargs={ 'extend':'neither' } )
        plt.title("Time = " + str(t2m.coords['time'].values[i])[:13])
        plt.savefig(F_png)
        plt.close()
os.system("ffmpeg -r 60 -f image2 -s 1920x1080 -i ./GRAPHICAL/animations/native_sequence/Python_Animation_03_frame_%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ./era5_native_t2m_2005_1st90days.mp4")
F_del = glob.glob('./GRAPHICAL/animations/native_sequence/*.png')
for F_delete in F_del: os.remove(F_delete)
for i in range(90*24):
    F_png = f"./GRAPHICAL/animations/regrid_sequence/Python_Animation_03_frame_{i:04}.png"
    if not os.path.exists(F_png):
        print('plotting: '+str(t2m.coords['time'].values[i])[:13])
        t2m.isel(time=i).t2m.plot( figsize = (12,6),
                                   vmin=190, vmax=340,
                                   cmap='coolwarm',    # Change the colormap back to 'bwr'
                                   cbar_kwargs={ 'extend':'neither' } )
        plt.title("Time = " + str(t2m.coords['time'].values[i])[:13])
        plt.savefig(F_png)
        plt.close()
os.system("ffmpeg -r 60 -f image2 -s 1920x1080 -i ../GRAPHICAL/animations/regrid_sequence/Python_Animation_03_frame_%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ./era5_regrid_t2m_2005_1st90days.mp4")
F_del = glob.glob('./GRAPHICAL/animations/regrid_sequence/*.png')
for F_delete in F_del: os.remove(F_delete)
