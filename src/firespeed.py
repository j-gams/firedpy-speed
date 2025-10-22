import math
import numpy as np
from shapely import LineString
import geopandas as gpd
import pyproj


def auto_density_func(density_param, shape, ii):
    ### thinking on this is we want the density per-square to remain similar...?
    ### so density to go with the sqrt of n_points
    ### arbitrarily set pt density to 1...
    n_pts = len(shape)
    max_bins = int(math.sqrt(n_pts)/density_param)
    ## print("bin density verification:", len(shape), "pts, ", len(shape)/(max_bins * max_bins), "LD", ii)

    return max_bins

def auto_simplify_func(bbox_lists):
    aggregate_bbox= (min([bbox_lists[ii][0] for ii in range(len(bbox_lists))]), min([bbox_lists[ii][1] for ii in range(len(bbox_lists))]),
                     max([bbox_lists[ii][2] for ii in range(len(bbox_lists))]), max([bbox_lists[ii][3] for ii in range(len(bbox_lists))]))

    ag_len_density = min(aggregate_bbox[2] - aggregate_bbox[0], aggregate_bbox[3] - aggregate_bbox[1])

    ### auto scale factor -- this can be hard coded?
    asf = 0.0001
    print("computed aggregate scale for simplification:", ag_len_density, ag_len_density * asf)

    return ag_len_density * asf

def firespeed_id_wrapper(fire_gdf, test_subset=float("inf")):
    unique_fireids = np.unique(fire_gdf["id"])
    for i in range(len(unique_fireids)):
        ### TODO -- need to only operate on two most recent days for each fire...?
        ### unsure about structure of the dataframe here
        fireid_gdf = fire_gdf[fire_gdf['id'] == unique_fireids[i]]
        


def computefirespeed(fire_gdf, test_subset=float("inf")):
    transformer = pyproj.Transformer.from_crs(fire_gdf.crs, "EPSG:4326", always_xy=True)
    geod = pyproj.Geod(ellps="WGS84")
    orig_x = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    orig_y = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    dest_x = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    dest_y = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    result_max_dist = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    result_speed = [np.nan for iii in range((min(fire_gdf.shape[0], 1)))]
    if len(fire_gdf) == 0:
        return orig_x, orig_y, dest_x, dest_y, result_max_dist, result_speed
    prev_step = [fire_gdf.iloc[0]["geometry"].geoms[ii].simplify(0.05).exterior.coords for ii in range(len(fire_gdf.iloc[0]["geometry"].geoms))]
    ### need to deal with resampling or something here in future iteration of code
    ### iterate over time steps
    for i in range(1, min(fire_gdf.shape[0], test_subset)):
        ### setup for overlap and spot checks...
        inter_matrix = np.zeros((len(fire_gdf.iloc[i-1]["geometry"].geoms), len(fire_gdf.iloc[i]["geometry"].geoms)))
        for ii in range(inter_matrix.shape[0]):
            for jj in range(inter_matrix.shape[1]):
                inter_matrix[ii, jj] = fire_gdf.iloc[i-1]["geometry"].geoms[ii].intersects(fire_gdf.iloc[i]["geometry"].geoms[jj])
        ### bug check
        if inter_matrix.shape[0] == 1 and inter_matrix.shape[1] == 1 and inter_matrix[0, 0] == False:
            print("Error: No overlap at", i)

        ### return maximum_distance, max_dist_origin, max_dist_destination
        max_fire_dist, max_origin, max_destination, prev_step = compute_max_vector(fire_gdf["geometry"][i-1].geoms, fire_gdf["geometry"][i].geoms, prev_step, 
                                                                                   inter_matrix, buffer=200, maxbins=200, slop=2)
        orig_x.append(max_origin[0])    ### error is here
        orig_y.append(max_origin[1])
        dest_x.append(max_destination[0])
        dest_y.append(max_destination[1])

        ### compute meter distance
        lons, lats = transformer.transform([max_origin[0], max_destination[0]],
                                           [max_origin[1], max_destination[1]])
        dist = geod.line_length(lons, lats)
        result_max_dist.append(dist)

        ### assumption here is these are daily perims...
        ### so to compute spread in km/h, we do (dist in km) / 24
        result_speed.append((dist * 1000) / 24)

    
        
    return orig_x, orig_y, dest_x, dest_y, result_max_dist, result_speed

def compute_max_vector(perim_inner_geoms, perim_outer_geoms, inner_coords, inter_matrix, bufferpct, maxbins, slop, simplify_param, iter_i, debug=False):
    ### distance computation 2 -- binned
    outer_coords = []

    ### remove later
    #buffer = bufferpct

    root2 = round(math.sqrt(2), 5)

    verify_outer = False

    result_over_polys = []

    ### perim_inner/outer geoms to be just gpd.iloc[i/i+1].geoms
    ### broadly, for points vi, vj and polys Px Py,
    ### minimize distances between vi and vj with i fixed (closest point to fixed point)
    ### maximize distances between vi and vjmax with vjmax the closest point to vi (furthest travelled between 2 polys)
    ### minimize distances between Px and Py with x fixed (closest polygon to the one looked at)
    ### maximize distances between Px and Pymax (furthest point from closest point on closest polygon)
    for poly_outer in range(len(perim_outer_geoms)):
        ### computer poly_outer bounding box
        outer_bbox = perim_outer_geoms[poly_outer].exterior.bounds
        outer = perim_outer_geoms[poly_outer].simplify(simplify_param).exterior.coords
        outer_coords.append(outer)

        ### this checks if any previous perimeter is inside this one...
        ### if not, spotted
        spot_flag = not np.any(inter_matrix[:, poly_outer])
        polyids = []
        if spot_flag:
            ### in the R code, this is resampling to get more points per meter to get accurate estimate...?
            ### need to compare with all polygons
            polyids= [ii for ii in range(len(perim_inner_geoms))]
        ### if it hasn't spotted, compare to all perims that fall inside:
        else:
            for ii in range(len(perim_inner_geoms)):
                if inter_matrix[ii, poly_outer]:
                    polyids.append(ii)

        
        ### now, keep track of all comparisons for this polygon:
        ### [distance, initial_coord, end_coord, pct_in_bounds, initial_poly, end_poly]
        poly_min = [float("inf"), None, None, None, None, None]

        verify_inner = False

        for poly_inner in polyids:
            if debug:
                print("working on poly combination", poly_outer, poly_inner)
            ### inner union:
            union_poly = shapely.unary_union([perim_outer_geoms[poly_outer], perim_inner_geoms[poly_inner]])
            ### compute inner bounding box...
            inner_bbox = perim_inner_geoms[poly_inner].exterior.bounds
            ### compute combined bounding box
            bbox = (min(inner_bbox[0], outer_bbox[0]), min(inner_bbox[1], outer_bbox[1]), 
                    max(inner_bbox[2], outer_bbox[2]), max(inner_bbox[3], outer_bbox[3]))
            ### guess x, y bin sizes 
            bins_x = math.ceil((bbox[2] - bbox[0])/maxbins)
            bins_y = math.ceil((bbox[3] - bbox[1])/maxbins)
            ### bin resolution is bigger of these values since we want square bins
            bin_res = max(bins_x, bins_y)
            ### record this param somewhere...
            adparam = 0.1
            bin_res = math.ceil(max(bbox[2] - bbox[0], bbox[3] - bbox[1]) / min(auto_density_func(adparam, perim_outer_geoms[poly_outer].simplify(simplify_param).exterior.coords, iter_i),
                                                                        auto_density_func(adparam, perim_inner_geoms[poly_inner].exterior.coords, iter_i)))
            
            ### compute actual number of bins with bin resolution 
            grid_size = (math.ceil((bbox[2] - bbox[0])/bin_res), math.ceil((bbox[3] - bbox[1])/bin_res))
            if debug:
                print("- debug -", grid_size)
            grid_spatial = (grid_size[0] * bin_res, grid_size[1] * bin_res)
            grid_offset = ((grid_spatial[0] - (bbox[2] - bbox[0]))/2, (grid_spatial[1] - (bbox[3] - bbox[1]))/2)

            ### now make np array w/ coarser side...
            inner_bins = np.zeros(grid_size, dtype=object)
            inner_occu = np.zeros(grid_size, dtype=bool)
            outer_bins = np.zeros(grid_size, dtype=object)
            outer_occu = np.zeros(grid_size, dtype=bool)
            for i in range(grid_size[0]):
                for j in range(grid_size[1]):
                    inner_bins[i, j] = [] 
                    outer_bins[i, j] = []
            ### lower left of grid
            lower = (bbox[0] - grid_offset[0], bbox[1] - grid_offset[1])

            ### now, bin inner layer
            for i in range(len(inner_coords[poly_inner])):
                ### bin inner layer 
                inner_ids = (int((inner_coords[poly_inner][i][0]-lower[0])//bin_res), 
                           (int(inner_coords[poly_inner][i][1]-lower[1])//bin_res))
                inner_bins[inner_ids[0], inner_ids[1]].append(i)
                inner_occu[inner_ids[0], inner_ids[1]] = True
            ### now, bin outer layer
            for i in range(len(outer)):
                ### bin outer layer
                outer_ids = (int((outer[i][0]-lower[0])//bin_res), (int(outer[i][1]-lower[1])//bin_res))
                outer_bins[outer_ids[0], outer_ids[1]].append(i)
                outer_occu[outer_ids[0], outer_ids[1]] = True
            
            ### compute the list of ids of occupied outer bins...
            outer_where = np.argwhere(outer_occu == True)

            all_bins_min = [float("-inf"), None, None, None]

            if debug:
                #print("- debug -", inner_bins, outer_bins)
                #print("- debug -", inner_occu, outer_occu)
                print("- debug: occupied outer bins", len(outer_where))

            nothing_found_flag = True

            ### finally ... we can do binned nearest neighbors
            ### do root2 rings... focus on points in outer perim
            ### iterate over occupied bins in outer ring to narrow comparison
            for occu_loc in outer_where:
                ### TODO -- GENERAL CLEAN UP
                ### BINNED METHOD -- 
                ### start at the center bin and iteratively look further until we find all squares...
                ### ... within 2sqrt(2) + slop of the closest point we find

                ### ...because there are only so many bins... loose upper bound with this:
                ring_bound = max(grid_size)
                ring_sqs = None
                ring_offset = None

                ### BIG NOTE:
                ### NEED TO DIFFERENTIATE BETWEEN
                ### - the case where no points are found and the case where points are found that are not valid

                temp_pts_dict = {}

                ### need to reformulate this algorithm to the following:
                ### --- look for first ring F
                ### --- compare all points in [F, F+root2offset]
                ### --- if we find something good to go... we have found a legal sample
                ### --- if we don't find something it is time to loop back
                ### --- go to bounds...
                ### --- wrap it up
                pair_found = False
                failed_prev = False
                ring1 = 0
                ring0 = 0

                ### CONCEPT
                ### we check at the ring for legal samples
                ### if we find samples but they are all illegal, we need to mask them out and keep trying
                
                ### this is the loop to find legal samples
                #while no_lower_found:
                ### step 1: look for ring with default parameters
                ### this loop is making sure the ring doesn't go to inf
                while ring1 < ring_bound:
                    ### root2_dist = math.ceil(root2 * (ring + 1)) + params["slop"]
                    root2_bound = math.ceil(root2 * (ring1 + 1)) + slop
                    ### in this case, we have failed to find a legal sample in a previous step.
                    if failed_prev:
                        if debug:
                            print("setting up ring1 mask due to previous failed attempt..")
                        ### idea here is to create a boolean mask of the size of the new outer ring...
                        ### need to bound this too with the bounds of the grid as a whole
                        temp1_mask = np.zeros((min(occu_loc[0] + ring1 + 1, grid_size[0]) - max(occu_loc[0] - ring1, 0), 
                                                min(occu_loc[1] + ring1 + 1, grid_size[1]) - max(occu_loc[1] - ring1, 0)), dtype=np.bool)
                        
                        ### set the inner ring (previous attempt area) to 1, indicating it should be masked out
                        temp1_mask[min(occu_loc[0], ring1) - min(occu_loc[0], ring0) : min(grid_size[0] - occu_loc[0], ring1) + min(grid_size[0] - occu_loc[0], ring0),
                                    min(occu_loc[0], ring1) - min(occu_loc[0], ring0) : min(grid_size[0] - occu_loc[0], ring1) + min(grid_size[0] - occu_loc[0], ring0)] = 1
                        
                        ### get coordinates of in-bounds occupied bins where the mask ring is 0 (do not mask out)
                        temp1 = np.where((inner_occu[max(occu_loc[0] - ring1, 0): min(occu_loc[0] + ring1 + 1, grid_size[0]), 
                                                        max(occu_loc[1] - ring1, 0): min(occu_loc[1] + ring1 + 1, grid_size[1])]) & (temp1_mask == 0))
                        
                        ### with this we know what bins to iterate over
                    ### if we have not yet failed to find legal samples (this is the first step..)
                    else:
                        ### just do the outer bounding ring
                        temp1 = np.where(inner_occu[max(occu_loc[0] - ring1, 0): min(occu_loc[0] + ring1 + 1, grid_size[0]), 
                                                    max(occu_loc[1] - ring1, 0): min(occu_loc[1] + ring1 + 1, grid_size[1])] == True)
                        
                    ### temp1 is [x coord list, y coord list]
                    ### if we have found a bin containing a point, then we can check for shortest with legality..
                    if debug:
                        if len(temp1[0]) > 0:
                            pass
                            #print("temp1 length:", len(temp1[0]), temp1, )
                    if len(temp1[0]) > 0:
                        if debug:
                            print("   - found", len(temp1[0]), "bins in ring1", len(outer_bins[occu_loc[0], occu_loc[1]]), ring1)
                        ### this probably needs a review ... general idea is that if this is a spot fire, 
                        ### we do NOT need to verify if the path to the point lies within bounds...

                        ### iterate over all points in the outer perim center grid sq.
                        for k in range(len(outer_bins[occu_loc[0], occu_loc[1]])):
                            ### remember what the point of this is..
                            ### for this specific outer occu_loc, set up a dict to store the comparison fields
                            if k not in temp_pts_dict:
                                ### idea of this is to capture
                                ### [distance, start_coords, end_coords, pct_within]
                                temp_pts_dict[k] = [float("inf"), None, None, None, False]
                            ### iterate over all located inner occupied bins
                            for i in range(len(temp1[0])):
                                ### extract contents of occupied bin at index i
                                temp_box = inner_bins[temp1[0][i] + max(occu_loc[0]-ring1, 0),
                                                        temp1[1][i] + max(occu_loc[1]-ring1, 0)]
                                ### iterate over pionts in occupied bin
                                for l in range(len(temp_box)):
                                    ### new strategy
                                    direct_line = LineString([(inner_coords[poly_inner][temp_box[l]]),
                                                                outer[outer_bins[occu_loc[0], occu_loc[1]][k]]])
                                    
                                    ### compute direct line length
                                    direct_line_length = direct_line.length
                                    ### is this true: SAVE THIS DIRECT LINE LENGTH FOR EACH POINT
                                    ### is this true: IF WE CAN't FIND SOMETHING WITHIN BUFFER OF THIS LENGTH MEETING 75% WITHIN, CANCEL!
                                    
                                    ### to optimize, compare the best-case scenario for this point (that it is legal) 
                                    ### ... with the nearest inner point yet found for this outer point.
                                    ### ... if it is less, we need to check that it is legal
                                    if direct_line_length < temp_pts_dict[k][0]:
                                        ### check legality
                                        ### if the line length is nonzero
                                        if direct_line_length > 0:
                                            ### TODO -- need to verify it lies within
                                            if spot_flag:
                                                ### TODO: NEED TO VERIFY IT ISN'T TOO FAR??
                                                ### Danielle / Nicole: What is the new condition for checking if vectors to spots are OK?
                                                temp_pts_dict[k][0] = direct_line_length
                                                temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                                temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                                temp_pts_dict[k][3] = 0
                                                temp_pts_dict[k][4] = True
                                                ring_bound = root2_bound
                                                pair_found = True
                                            else:
                                                ### intersection...? what is this doing...
                                                direct_inter = direct_line.intersection(union_poly)
                                                ### compute the intersecting length with this... if there is no intersect, 0
                                                intersection_length = direct_inter.length if not direct_inter.is_empty else 0.0
                                                ### check whether the intersecting length meets the threshold 
                                                if intersection_length / direct_line_length >= bufferpct:
                                                    ### we have caught a legal sample, so we now know the upper bound is r2d
                                                    ### ... which we already computed. Save time by finishing the loop
                                                    ### ... and continuing to work on temp2 below
                                                    pair_found = True
                                                    temp_pts_dict[k][0] = direct_line_length
                                                    temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                                    temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                                    temp_pts_dict[k][3] = intersection_length / direct_line_length
                                                    temp_pts_dict[k][4] = True
                                                    ring_bound = root2_bound
                                        ### if the line length is 0
                                        else:
                                            ### we found a line with length 0... this is necessarily a match
                                            ### do we need to ensure this is within some very small epsilon?
                                            pair_found = True
                                            temp_pts_dict[k][0] = direct_line_length
                                            temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                            temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                            temp_pts_dict[k][3] = 1
                                            temp_pts_dict[k][4] = True
                                            ring_bound = root2_bound
                        if debug:
                            print("   - done with ring", ring1, pair_found)
                        ### IF WE DID NOT FIND ANY LEGAL COMBINATIONS IN THE FIRST RING:
                        ### this is the case where we found points but they are not legal (oob)
                        ### we need to set a mask so that we don't redo the work of checking in-bounds
                        ### ...then increment the ring size and skip to next iteration
                        if not pair_found:
                            failed_prev = True
                            ring0 = ring1
                            ring1 += 1
                            continue
                        ### TODO -- is there a further optimization here to increase the size of the mask to exclude temp1?
                        ### this is the same deal but for the outer ring at this step so we don't have to do 2 steps.
                        temp2_mask = np.zeros((min(occu_loc[0] + root2_bound + 1, grid_size[0]) - max(occu_loc[0] - root2_bound, 0), 
                                                min(occu_loc[1] + root2_bound + 1, grid_size[1]) - max(occu_loc[1] - root2_bound, 0)), dtype=np.bool)
                        ### this time, mask out ring1 as well to avoid duplicating computation
                        temp2_mask[min(occu_loc[0], root2_bound) - min(occu_loc[0], ring1) : min(grid_size[0] - occu_loc[0], root2_bound) + min(grid_size[0] - occu_loc[0], ring1),
                                    min(occu_loc[0], root2_bound) - min(occu_loc[0], ring1) : min(grid_size[0] - occu_loc[0], root2_bound) + min(grid_size[0] - occu_loc[0], ring1)] = 1
                        temp2 = np.where((inner_occu[max(occu_loc[0] - root2_bound, 0): min(occu_loc[0] + root2_bound + 1, grid_size[0]), 
                                                        max(occu_loc[1] - root2_bound, 0): min(occu_loc[1] + root2_bound + 1, grid_size[1])]) & (temp2_mask == 0))
                        ### set exit condition here to be sure because we know we have already found a pair in ring1
                        ring1 = root2_bound
                        if len(temp2[0]) > 0:
                            ### if we found a legal pair above we can jump straght back in with temp2
                            for k in range(len(outer_bins[occu_loc[0], occu_loc[1]])):
                                    ### iterate over all located inner occupied bins
                                for i in range(len(temp2[0])):
                                    ### extract contents of occupied bin at index i
                                    temp_box = inner_bins[temp2[0][i] + max(occu_loc[0]-root2_bound, 0),
                                                            temp2[1][i] + max(occu_loc[1]-root2_bound, 0)]
                                    ### iterate over pionts in occupied bin
                                    for l in range(len(temp_box)):
                                        ### new strategy
                                        direct_line = LineString([(inner_coords[poly_inner][temp_box[l]]),
                                                                    outer[outer_bins[occu_loc[0], occu_loc[1]][k]]])
                                        direct_line_length = direct_line.length
                                        if direct_line_length < temp_pts_dict[k][0]:
                                            if direct_line_length > 0:
                                                if spot_flag:
                                                    ### TODO: NEED TO VERIFY IT ISN'T TOO FAR? what is the spot fire condition
                                                    pair_found = True
                                                    temp_pts_dict[k][0] = direct_line_length
                                                    temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                                    temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                                    temp_pts_dict[k][3] = 0
                                                    temp_pts_dict[k][4] = True
                                                else:
                                                    direct_inter = direct_line.intersection(union_poly)
                                                    ### compute the intersecting length with this... if there is no intersect, 0
                                                    intersection_length = direct_inter.length if not direct_inter.is_empty else 0.0
                                                    ### check whether the intersecting length meets the threshold 
                                                    if intersection_length / direct_line_length >= bufferpct:
                                                        ### we have caught a legal sample, so we now know the upper bound is r2d
                                                        ### ... which we already computed. Save time by finishing the loop
                                                        ### ... and continuing to work on temp2 below
                                                        pair_found = True
                                                        temp_pts_dict[k][0] = direct_line_length
                                                        temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                                        temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                                        temp_pts_dict[k][3] = intersection_length / direct_line_length
                                                        temp_pts_dict[k][4] = True
                                            else:
                                                pair_found = True
                                                temp_pts_dict[k][0] = direct_line_length
                                                temp_pts_dict[k][1] = inner_coords[poly_inner][temp_box[l]]
                                                temp_pts_dict[k][2] = outer[outer_bins[occu_loc[0], occu_loc[1]][k]]
                                                temp_pts_dict[k][3] = 1
                                                temp_pts_dict[k][4] = True
                            ### in this case we have found legal elements in the first ring, extracted an outer ring ... 
                            ### ... with the root2 formula, so now we have at least some nearest neighbor.
                            ### ... so we can set ring to root2_dist and end the iteration
                            ring1 = root2_bound

                    ### This is the case where we do not find anything at all, so we can safely just increment ring1...
                    ### ... and not deal with the masking
                    else:
                        ring1 += 1
                        continue

                ### we found a pair...
                if pair_found:
                    nothing_found_flag = False
                ### iteration is complete, up to ring_bound, and no pair has been found (uh oh)
                ### go to next step before we error
                else:
                    if debug:
                        print("- debug - no rings found on this occu_loc", ring0, ring1, root2_bound)
                    continue
                
                
                bin_min = [float("-inf"), None, None, None]
                ### get the shortest pair here..
                ### set bin_min to the shortest pair found for this entire outer bin
                for k in range(len(outer_bins[occu_loc[0], occu_loc[1]])):
                    ### break ties with whichever has a higher percentage within the perim..?
                    if temp_pts_dict[k][4] and temp_pts_dict[k][0] > bin_min[0] or (temp_pts_dict[k][0] == bin_min[0] and temp_pts_dict[k][3] > bin_min[3]):
                        bin_min = temp_pts_dict[k][:4]
                ### now we can compare the longest shortest-pair combo we found from this bin to the ones
                ### found in other bins...
                if bin_min[3] is None and debug:
                    print(bin_min, pair_found, ring1, len(outer_bins[occu_loc[0], occu_loc[1]]))
                    print(temp_pts_dict)
                if bin_min[0] > all_bins_min[0] or (bin_min[0] == bin_min[0] and bin_min[3] > all_bins_min[3]):
                    all_bins_min = bin_min
            if nothing_found_flag:
                continue
            ### at this point we have the furthest nearest-neighbor for this polygon pair
            if all_bins_min[0] < poly_min[0] or (all_bins_min[0] == poly_min[0] and all_bins_min[3] > poly_min[3]):
                poly_min = all_bins_min + [poly_inner, poly_outer]
                ### set this to true to carry through to the end that yes, we have found some legal coord pair somewhere, 
                ### ... so we will not need to return NaNs or whatever
                verify_inner = True
        ### aggregate over all outer polygons
        if verify_inner:
            result_over_polys.append(poly_min)
            ### set this to true to carry through to the end that yes, we have found some legal coord pair somewhere, 
            ### ... so we will not need to return NaNs or whatever
            verify_outer = True
        else:
            continue
    if verify_outer:
        ### gather only the lenghts
        rop_np = np.array([result_over_polys[ii][0] for ii in range(len(result_over_polys))])
        max_dist_pair = np.argmax(rop_np)
        maximum_distance = result_over_polys[max_dist_pair][0]
        max_dist_origin = result_over_polys[max_dist_pair][1]
        max_dist_destination = result_over_polys[max_dist_pair][2]
        max_pct_inbounds = result_over_polys[max_dist_pair][3]
        max_poly_origin = result_over_polys[max_dist_pair][4]
        max_poly_destination = result_over_polys[max_dist_pair][5]
    else:
        ### debug plots were here
        max_loc = None
        maximum_distance = np.nan
        max_dist_origin = None
        max_dist_destination = None
        max_poly_origin = None
        max_poly_destination = None
    return maximum_distance, max_dist_origin, max_dist_destination, outer_coords, max_pct_inbounds, max_poly_origin, max_poly_destination

### fire distance computations
### direct (as crow flies) distance
def compute_dist(a, b):
    return math.sqrt(((a[0] - b[0]) ** 2) + ((a[1] - b[1]) ** 2))