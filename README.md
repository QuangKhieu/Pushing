# Multi-Robot Cooperative Object Pushing

[link_demo:](https://youtu.be/hgaxgAIuazg)

This project develops a distributed control strategy for a multi-robot formation to autonomously approach, cooperatively push a heavy object, and automatically disengage after completing the task without direct communication.

## Behavior Control Mechanism
- Each robot autonomously understands its role in the task.
- The formation finds optimal pushing positions to overcome frictional forces, move in the correct direction, and optimize speed.

## Adaptive Task Mechanism
- Uses finite state machines to enable the formation to flexibly switch behaviors.

## Simulation and Experiments
- **Platform**: Simulated in the Vrep physical environment.
- **Robot Design**: 4-wheel omni-directional robots.
- **Experiments**: Tested the control strategy in various scenarios with different environmental and physical parameters.

## Results
- The strategy meets the criteria for distribution, adaptability, and flexibility.
- Demonstrates potential for real-world application.
