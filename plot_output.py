import pickle, lzma
import plot.plots as plotting

filename = '/Users/leecarlaw/Desktop/test/sharppy.xz'
with lzma.open(filename, 'rb') as f: data = pickle.load(f)
plotting.create_map_output(data[0], data[0]['lons'], data[0]['lats'])