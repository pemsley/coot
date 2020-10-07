
#include "matrix-utils.hh"

clipper::Mat33<double> glm_to_mat33(const glm::mat4 &quat_mat) {

   clipper::Mat33<double> m(quat_mat[0][0], quat_mat[0][1], quat_mat[0][2],
                            quat_mat[1][0], quat_mat[1][1], quat_mat[1][2],
                            quat_mat[2][0], quat_mat[2][1], quat_mat[2][2]);
   return m;

}

glm::quat coot_quaternion_to_glm(const coot::util::quaternion &q) {

   return glm::quat(q.q0, q.q1, q.q2, q.q3);

}

coot::util::quaternion glm_to_coot_quaternion(const glm::quat &q) {

   coot::util::quaternion qq(q[3], q[0], q[1], q[2]);
   return qq;
}
