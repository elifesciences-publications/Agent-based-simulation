# Agent-based-simulation
Agent based simulation of constant and active disassembly processes to simulate length regulation in Chlamydomonas.

Each program takes parameters which the user can change-: time over which the process is simulated in minutes(n) Diffusion constant(D), total tubulin in micron(T), gamma in the model(gammma), value of gamma*k*M(gammakM), constant disassembly constant term(da), velovity of ballistic motion(v), steady state length(Lss), time at which flagella is severed(in minutes), time step(dt), (replenishment time for motors and tubulins(t_m and t_t respectively);only  in TM_shared_depoly)

Output of each program is a plot of length of each flagella vs time.

Files include-:

TM_shared_depoly.m-: Agent based simulation of active disassembly process with the option of including feedback and severing.

two_flagella_motor_diff_tubulins_diff.m-: Agent based simulation of constant disassembly process with separate motors and tubulin pools with the option of including severing.

two_flagella_motor_diff_tubulins_shared.m-: Agent based simulation of constant disassembly process with separate motors pools but shared tubulin pool with the option of including severing.

two_flagella_motor_shared_tubulins_diff.m-: Agent based simulation of constant disassembly process with shared motors pool but separate tubulin pools with the option of including severing.

two_flagella_motor_shared_tubulins_shared.m-: Agent based simulation of constant disassembly process with shared motors and tubulin pool with the option of including severing.

### This code is associated with the paper from Fai et al., "Length regulation of multiple flagella that self-assemble from a shared pool of components". eLife, 2019. http://dx.doi.org/10.7554/eLife.42599

