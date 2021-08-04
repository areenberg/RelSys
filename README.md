# RelSys
RelSys is a tool for evaluating a multiclass system of M/M/c/c queues that are connected *only* through customer relocations (i.e. customers that are transferred to an alternative queue instead of being rejected from the system).
The tool is written in C++ and currently employs two different approaches for evaluating the system. The first is a heuristic mathematical model, which is based on a continuous-time Markov chain, and the second is a discrete-event simulation. The setup of the system is the same in both cases, but each approach comes with different methods for evaluating the system and viewing the results.

# Table of contents

1. How does it work
2. Getting started
3. How to cite
4. Licence

# How does it work

Consider a number of parallel M/M/c/c queues. That is, queues where customers arrive according to a Poisson process and have exponentially distributed service-time. In the common M/M/c/c queue, which is often denoted an Erlang loss or Erlang-B system, a customer is rejected (and lost from the system) if there are no idle servers upon arrival to the queue. However, in the RelSys-modeling tool, we allow customers to be transferred to one of the other queues with known probability. If there is an idle server in the alternative queue, the customer is served with an exponentially distributed time with the same rate-parameter as the customer would have had in the original queue. Thus, the system is *multiclass*. 


<img src="https://github.com/areenberg/RelSys/blob/development/images/example_system.jpeg" width="435" height="90">

# Getting started

This is how to get started


# How to cite

Coming soon.

# License

Copyright 2021 Anders Reenberg Andersen.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
