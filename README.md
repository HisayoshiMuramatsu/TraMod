# TraMod - Position Trajectory Modification Algorithm for Position, Velocity, and Acceleration Constraints

This is a real-time position trajectory modification algorithm in C++. It modifies a given position trajectory with position, velocity, and acceleration constraints. This algorithm was used by a mobile quad-arm robot at each joint not to violate its practical angle, angular velocity, and angular acceleration constraints and is introduced in the paper[1].

[1] Hisayoshi Muramatsu, Keigo Kitagawa, Jun Watanabe, Yuika Yoshimoto, and Ryohei Hisashiki, “A Mobile Quad-Arm Robot ARMS: Wheeled-Legged Tripedal Locomotion and Loco-Manipulation,” Journal of Robotics and Mechatronics, vol. 37, no. 2, pp. 489-499, Apr. 2025.
https://www.fujipress.jp/jrm/rb/robot003700020489/

## Example

There is a program that tests the position trajectory modification algorithm. Note that this algorithm needs Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page).

## Licence

[MIT License](https://github.com/HisayoshiMuramatsu/PASF/blob/master/LICENSE) © Hisayoshi Muramatsu
