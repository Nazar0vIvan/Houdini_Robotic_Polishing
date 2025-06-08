// VEX implementation for Houdini
function int[] findAvailableIndices(int num_points) {
    int indices[];
    for (int i = 0; i < num_points; i++) {
        append(indices, i);
    }
    return indices;
}

function float planeDistance(vector p; float A, B, C, D) {
    return A*p.x + B*p.y + C*p.z + D;
}

function int[][] findPointPairs(vector curves[][]; float A, B, C, D) {
    int pairs[][];
    int num_curves = len(curves);
    if (num_curves < 2) return pairs;
    
    int num_points = len(curves[0]);
    int available_indices[] = findAvailableIndices(num_points);
    
    for (int j = 0; j < num_curves - 1; j++) {
        int i = 0;
        while (i < len(available_indices)) {
            int idx = available_indices[i];
            vector p_below = curves[j][idx];
            vector p_above = curves[j+1][idx];
            
            float dist_below = planeDistance(p_below, A, B, C, D);
            float dist_above = planeDistance(p_above, A, B, C, D);
            
            if (dist_below < 0 && dist_above > 0) {
                // Found a valid pair
                int pair[] = array(idx, j, j+1); // Stores (point_idx, curve_below, curve_above)
                append(pairs, pair);
                
                // Remove this index from future consideration
                removeindex(available_indices, i);
                
                // Early exit if all points paired
                if (len(pairs) == num_points) {
                    return pairs;
                }
            } else {
                i++;
            }
        }
    }
    return pairs;
}
