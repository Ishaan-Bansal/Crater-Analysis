import ray_tracing_curvature as rtc
import ray_tracing_X_Slices as rtx
import ray_tracing_Z_Slices as rtz

# # Write the relative path to the folder you want to load
# path = "Crater_STL_Files"
# cwd = os.getcwd()
# os.chdir(path)
# for filename in os.listdir():
#     if filename[-4:] != ".stl":
#         continue
#     locations = z_slice(filename, 20)
#     x, y = locations[:, 0], locations[:, 1]
#     plt.plot(x, y, 'k.', linewidth=2)
#     plt.axis('equal')
#     # plt.show()
#     filename = filename[:-4]
#     os.chdir(cwd)
#     plt.savefig("Crater_Slices/" + filename + "_Z-Slice" + ".png")
#     os.chdir(path)
#     plt.close()