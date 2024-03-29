# TraMod - Position Trajectory Modification Algorithm for Position, Velocity, and Acceleration Constraints

This is a real-time position trajectory modification algorithm in C++. It modifies a given position trajectory with position, velocity, and acceleration constraints. This algorithm was used by a mobile quad-arm robot at each joint not to violate its practical angle, angular velocity, and angular acceleration constraints and is introduced in the paper[1].

[1] Hisayoshi Muramatsu, Keigo Kitagawa, Jun Watanabe, Yuika Yoshimoto, and Ryohei Hisashiki, “A Mobile Quad-Arm Robot ARMS: Wheeled-Legged Tripedal Locomotion and Quad-Arm Loco-Manipulation,” arXiv, arXiv:2305.01406, May 2023. https://arxiv.org/abs/2305.01406

## Example

There is a program that tests the position trajectory modification algorithm. Note that this algorithm needs Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page).

## Licence

[MIT License](https://github.com/HisayoshiMuramatsu/PASF/blob/master/LICENSE) © Hisayoshi Muramatsu
