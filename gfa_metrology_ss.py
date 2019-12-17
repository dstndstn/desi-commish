from astrometry.util.fits import fits_table
import numpy as np

# Values copy-n-pasted from GFA Metrology spreadsheets.
# Specifically, "spreadsheet_*" are from the excel file named like
#  GFA10_Lateral_Measurements.xlsx
# in the "GFA images" tab.
# The "gif_ss_*_[12]" is from the "GIFS" tab.

# And we also unzip the images in a local directory --
# https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=4832;filename=GFA%2310_GIF-Metrology-images.zip;version=2

# DESI-4832-v2
spreadsheet_10 = '''		IMG_40.fits	Background reference	47.4		
38.5	58	IMG_41.fits	Quadrant 3top right	47.6	2112	537
50	58	IMG_42.fits	Quadrant 3top left	47.5	1346	534
56.5	58	IMG_43.fits	Quadrant 4 upper right	47.5	849	533
67.5	58	IMG_44.fits	Quadrant 4 upper left	47.4	116	530
45.5	63	IMG_45.fits	Quadrant 3 lower middle	47.5	1645	869
60.75	63	IMG_46.fits	Quadrant 4 lower middle	47.4	565	866
60.75	67	IMG_47.fits	Quadrant 1 upper middle	47.4	564	1180
45.5	67	IMG_48.fits	Quadrant 2 Upper Middle	47.5	1644	1184
38.5	71	IMG_49.fits	Quadrant 2 lower right	47.6	2110	1451
50	71	IMG_50.fits	Quadrant 2 lower left	47.6	1344	1449
56.5	71	IMG_51.fits	Quadrant 1 lower right	47.5	846	1448
67.5	71	IMG_52.fits	Quadrant 1 lower left	47.5	114	1445'''

# 	First Measurement			Second Measurement		
# Pinhole	X	Y	Z	X	Y	Z
gif_ss_10_1 = '''1	38.302	27.098	14.918			14.928
2	39.894	27.263	14.91			14.9
3	38.184	28.301	14.92			14.915
4	39.039	27.78	14.91			14.918'''
gif_ss_10_2 = '''1	78.236	75.348	15.035			15.028
2	79.655	74.631	15.01			15.02
3	78.769	76.417	15.025			15.035
4	79.217	75.526	15.04			15.03'''

# DESI-4823-v2
spreadsheet_6 = '''38.5	58	IMG_24.fits	Quadrant 3top right	46.2	2111	543
50	58	IMG_25.fit	Quadrant 3top left	46.5	1345	544
56.5	58	IMG_26.fit	Quadrant 4 upper right	46.2	848	545
67.5	58	IMG_27_1.fits	Quadrant 4 upper left	46.4	115	545
45.5	63	IMG_28.fits	Quadrant 3 lower middle	46.5	1644	878
60.75	63	IMG_29_1.fits	Quadrant 4 lower middle	46.3	565	878
60.75	67	IMG_30_1.fits	Quadrant 1 upper middle	47.1	565	1193
45.5	67	IMG_31_1.fits	Quadrant 2 Upper Middle	47.1	1645	1192
38.5	71	IMG_32_1.fits	Quadrant 2 lower right	46.7	2112	1458
50	71	IMG_33_1.fits	Quadrant 2 lower left	46.8	1346	1459
56.5	71	IMG_34.fits	Quadrant 1 lower right	47.3	849	1459
67.5	71	IMG_35.fits	Quadrant 1 lower left	47	116	1459
		IMG_37.fits	Background reference	48		'''

gif_ss_6_1 = '''1	39.561	28.842	14.9			14.9
2	38.018	28.495	14.9			14.898
3	39.833	27.663	14.89			14.9
4	38.926	28.073	14.9			14.895'''
gif_ss_6_2 = '''1	80.076	75.375	15.06			15.068
2	78.874	76.428	15.08			15.06
3	79.276	74.475	15.075			15.065
4	79.077	75.454	15.075			15.08'''

# DESI-4821-v2
spreadsheet_2 = '''38.5	58	IMG_51.fits	Quadrant 3top right	50.9	2112	547
52.5	58	IMG_54.fits	Quadrant 3top left	50.8	1180	544
54	58	IMG_56.fits	Quadrant 4 upper right	51	1016	544
67.5	58	IMG_58.fits	Quadrant 4 upper left	50.4	116	540
45.5	63	IMG_60.fits	Quadrant 3 lower middle	51	1645	879
60.75	63	IMG_61.fits	Quadrant 4 lower middle	50.6	564	875
60.75	67	IMG_62.fits	Quadrant 1 upper middle	51.3	564	1190
45.5	67	IMG_63.fits	Qandrant 2 Upper Middle	50.8	1644	1194
38.5	71	IMG_64.fits	Quadrant 2 lower right	51.2	2110	1462
52.5	71	IMG_65.fits	Quadrant 2 lower left	51.2	1178	1459
54	71	IMG_66.fits	Quadrant 1 lower right	50.6	1013	1548
67.5	71	IMG_67.fits	Quadrant 1 lower left	50.7	114	1455'''

gif_ss_2_1 = '''1	39.756	27.294	14.9	39.755	27.297
2	39.563	28.89	14.89	39.563	28.891
3	38.557	27.15	14.91	38.557	27.152
4	39.068	28.017	14.88	39.067	28.017'''
gif_ss_2_2 = '''1	79.573	74.311	15.015	79.574	74.311
2	79.942	75.872	15.02	79.94	75.873
3	78.41	74.584	15.025	78.411	74.583
4	79.17	75.226	15.02	79.174	75.225'''

# DESI-4836-v1
spreadsheet_8 = '''		IMG_22.fits	Background reference			
38.5	58	IMG_23.fits	Quadrant 3top right	50.1	2113	548
50	58	IMG_24.fits	Quadrant 3top left	49.7	1348	545
56.5	58	IMG_26.fits	Quadrant 4 upper right	49.7	850	544
67.5	58	IMG_27.fits	Quadrant 4 upper left	49.6	117	541
45.5	63	IMG_28.fits	Quadrant 3 lower middle	50.5	1646	880
60.75	63	IMG_29.fits	Quadrant 4 lower middle	49.9	566	876
60.75	67	IMG_30.fits	Quadrant 1 upper middle	50.1	565	1191
45.5	67	IMG_31.fits	Quadrant 2 Upper Middle	50.2	1645	1194
38.5	71	IMG_32.fits	Quadrant 2 lower right	49.9	2111	1463
50	71	IMG_33.fits	Quadrant 2 lower left	50.2	1345	1460
56.5	71	IMG_34.fits	Quadrant 1 lower right	50.3	848	1459
67.5	71	IMG_35.fits	Quadrant 1 lower left	49.9	115	1455'''

gif_ss_8_1 = '''1	38.865	28.995	14.945			14.945
2	37.957	27.677	14.952			14.932
3	39.849	28.311	14.94			14.952
4	38.9	27.991	14.935			14.93'''
gif_ss_8_2 = '''1	78.429	74.214	15.048			15.062
2	79.715	75.162	15.03			15.042
3	77.723	75.188	15.07			15.08
4	78.714	75.176	15.058			15.055'''

# DESI-4780-v2
spreadsheet_1 = '''38.5	58	IMG_18.fits	Quadrant 3top right	45	2118	572
52.5	58	IMG_20.fits	Quadrant 3top left	47	1186	573
54	58	IMG_23.fits	Quadrant 4 upper right	46	1021	574
67.5	58	IMG_28.fits	Quadrant 4 upper left	46	122	574
45.5	63	IMG_30.fits	Quadrant 3 lower middle	46	1651	906
60.75	63	IMG_32.fits	Quadrant 4 lower middle	46	572	907
60.75	67	IMG_34.fits	Quadrant 1 upper middle	47	572	1222
45.5	67	IMG_36.fits	Qandrant 2 Upper Middle	46	1651	1220
38.5	71	IMG_40.fits	Quadrant 2 lower right	46	2119	1486
52.5	71	IMG_42.fits	Quadrant 2 lower left	47	1187	1487
54	71	IMG_44.fits	Quadrant 1 lower right	47	1023	1488
67.5	71	IMG_46.fits	Quadrant 1 lower left	46	123	1488'''

gif_ss_1_1 = '''1	39.557	26.794	14.97	39.546	26.794	14.97
2	39.335	28.386	14.97	39.326	28.384	14.97
3	38.374	26.627	14.97	38.368	26.627	14.97
4	38.862	27.503	14.97	38.854	27.502	14.97'''
gif_ss_1_2 = '''1	80.148	74.608	15.075	80.14	74.606	15.075
2	79.064	75.783	15.075	79.059	75.784	15.075
3	79.255	73.786	15.075	79.251	73.785	15.075
4	79.158	74.791	15.075	79.155	74.79	15.075'''

# DESI-4822-v3
spreadsheet_4 = '''38.5	58	IMG_30.fits	Quadrant 3top right	49.5	2125	555
52.5	58	IMG_31.fits	Quadrant 3top left	49.1	1193	559
54	58	IMG_32.fits	Quadrant 4 upper right	49.5	1029	559
67.5	58	IMG_33.fits	Quadrant 4 upper left	48.9	129	562
45.5	63	IMG_34.fits	Quadrant 3 lower middle	48.9	1660	890
60.75	63	IMG_35.fits	Quadrant 4 lower middle	49.6	580	894
60.75	67	IMG_36.fits	Quadrant 1 upper middle	49	582	1209
45.5	67	IMG_37.fits	Qandrant 2 Upper Middle	49.4	1661	1205
38.5	71	IMG_38.fits	Quadrant 2 lower right	49.6	2129	1469
52.5	71	IMG_39.fits	Quadrant 2 lower left	49.1	1196	1473
54	71	IMG_40.fits	Quadrant 1 lower right	49.5	1033	1474
67.5	71	IMG_41.fits	Quadrant 1 lower left	49.1	133	1477'''

gif_ss_4_1 = '''1	39.226	26.915	14.87
2	39.595	28.482	14.89
3	38.071	27.189	14.91
4	38.837	27.833	14.88'''
##### NOTE that pinhole 1 here is wrong -- it's a re-measurement of pinhole 3.
gif_ss_4_2 = '''1	79.562	73.973	15.03
2	78.628	75.75	15.03
3	79.558	73.972	15.03
4	79.101	74.867	15.03'''

# THIS is a patched version -- I filled in pinhole 1 with a fake measurement
# --
# Pinhole 1 (row 0) is wrong.
# v = np.zeros((3,2))
# cx = G.x[3]
# cy = G.y[3]
# v[:,0] = G.x[1:] - cx
# v[:,1] = G.y[1:] - cy
# U,S,V = np.linalg.svd(v)
# fx = 0.96
# fy = 0.72
# v1x = G.x[2] - V[1,0]*fx - V[0,0]*fy
# v1y = G.y[2] - V[1,1]*fx - V[0,1]*fy
# 
gif_ss_4_2 = '''1	80.075	75.055	15.03
2	78.628	75.75	15.03
3	79.558	73.972	15.03
4	79.101	74.867	15.03'''

def parse_gfa_ss(ss_text):
    lines = [row.split('\t')
             for row in ss_text.split('\n')
             if len(row.strip())>0]
    # Fields are:
    # mmx, mmy, filename, name, fgpa_t, spotx, spoty
    G = fits_table()
    #print('names:', [w[3] for w in lines])
    G.is_bg = np.array(['Background' in w[3] for w in lines])
    G.mmx = np.array([float(w[0]) if w[0] != '' else -1
                      for w in lines]).astype(np.float32)
    G.mmy = np.array([float(w[1]) if w[1] != '' else -1
                      for w in lines]).astype(np.float32)
    G.img_fn = np.array([w[2] for w in lines])
    G.orig_ssx = np.array([float(w[5]) if w[5] != '' else -1
                           for w in lines]).astype(np.int16)
    G.orig_ssy = np.array([float(w[6]) if w[6] != '' else -1
                           for w in lines]).astype(np.int16)
    return G

def parse_gif_ss(ss):
    G = fits_table()
    lines = ss.split('\n')
    lines = [line.split('\t') for line in lines]
    G.pinhole = np.array([int(w[0]) for w in lines])
    G.x = np.array([float(w[1]) for w in lines])
    G.y = np.array([float(w[2]) for w in lines])

    def dist_btw(i, j):
        return np.hypot(G.x[i] - G.x[j], G.y[i] - G.y[j])
    def check_dist(d, target):
        if np.abs(d - target) >= 0.02:
            print('Distance supposed to be', target, 'but is', d)
        assert(np.abs(d - target) < 0.02)
    assert(np.all(G.pinhole == np.arange(1,5)))
    d12 = dist_btw(0, 1)
    check_dist(d12, 1.6)
    d13 = dist_btw(0, 2)
    check_dist(d13, 1.2)
    d14 = dist_btw(0, 3)
    check_dist(d14, 1.0)
    d23 = dist_btw(1, 2)
    check_dist(d23, 2.0)
    d24 = dist_btw(1, 3)
    check_dist(d24, 1.0)
    d34 = dist_btw(2, 3)
    check_dist(d34, 1.0)

    return G

gfa_spreadsheets = {10: parse_gfa_ss(spreadsheet_10),
                    6:  parse_gfa_ss(spreadsheet_6),
                    2:  parse_gfa_ss(spreadsheet_2),
                    8:  parse_gfa_ss(spreadsheet_8),
                    1:  parse_gfa_ss(spreadsheet_1),
                    4:  parse_gfa_ss(spreadsheet_4)
                    }

gif_spreadsheets = {
    1 :(parse_gif_ss(gif_ss_1_1),  parse_gif_ss(gif_ss_1_2)),
    2 :(parse_gif_ss(gif_ss_2_1),  parse_gif_ss(gif_ss_2_2)),
    4 :(parse_gif_ss(gif_ss_4_1),  parse_gif_ss(gif_ss_4_2)),
    6 :(parse_gif_ss(gif_ss_6_1),  parse_gif_ss(gif_ss_6_2)),
    8 :(parse_gif_ss(gif_ss_8_1),  parse_gif_ss(gif_ss_8_2)),
    10:(parse_gif_ss(gif_ss_10_1), parse_gif_ss(gif_ss_10_2)),}

