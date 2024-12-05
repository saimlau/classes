#!/usr/bin/env python3

import numpy as np
import rclpy

from asl_tb3_lib.control import BaseHeadingController
from asl_tb3_msgs.msg import TurtleBotControl, TurtleBotState
from asl_tb3_lib.math_utils import wrap_angle


class HeadingController(BaseHeadingController):

    def __init__(self) -> None:
        super().__init__()

        # control gains
        self.kp = 2.0

    def compute_control_with_goal(self,
        state: TurtleBotState,
        goal: TurtleBotState,
    ) -> TurtleBotControl:
        control = TurtleBotControl()

        err = wrap_angle(goal.theta - state.theta)
        control.omega = self.kp * err

        return control


if __name__ == "__main__":
    rclpy.init()
    node = HeadingController()
    rclpy.spin(node)
    rclpy.shutdown()