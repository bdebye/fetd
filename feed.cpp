
#include "feed.h"
#include "solve.h"

const int SOURCE_GROUP_TAG = 286;

int feed_index;

int find_feed_index() {
    msh_segment feed_segment;
    // std::cout << mesh.segment.size() << std::endl;
    for(int i = 0; i < mesh.segment_num; i++) {
        // std::cout << mesh.segment[i].physical_tag << std::endl;
        if(mesh.segment[i].physical_tag == SOURCE_GROUP_TAG) {
            feed_segment = mesh.segment[i];
            break;
        }
    }
    for(int i = 0; i < edge_basis.size(); i++) {
        if((edge_basis[i][0] == feed_segment.node_number_list[0] &&
            edge_basis[i][1] == feed_segment.node_number_list[1]) ||
           (edge_basis[i][0] == feed_segment.node_number_list[1] &&
            edge_basis[i][1] == feed_segment.node_number_list[0])) {
               return i;
           }
    }
    return -1;
}


void set_feed_index() {
    feed_index = find_feed_index();
    if(feed_index < 0) {
        std::cout << "The feed is not properly configured.." << std::endl;
    }
}