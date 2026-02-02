#include <ros/ros.h>
#include <tf/transform_broadcaster.h>

int main(int argc, char** argv){
  ros::init(argc, argv, "tf_mytrans");
  ros::NodeHandle n;

  ros::Rate r(100);//100hz

  tf::TransformBroadcaster broadcaster;//创建一个广播类 broadcaster

  while(n.ok()){
    tf::Transform transform;
    tf::Quaternion q;
    transform.setOrigin(tf::Vector3(0,0,0));
    q.setRPY(0, 0,M_PI/2);
    transform.setRotation(q);
    broadcaster.sendTransform(tf::StampedTransform(transform, ros::Time(0), "world", "trans"));
    r.sleep();                                    

  }
}

