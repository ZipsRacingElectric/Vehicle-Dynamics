% vcu_yaw_table Exports a yaw lookup table to program into the VCU.
%
%   vcu_yaw_table (A) - Prints the yaw lookup table defined by the matrix A.
%   Where:
%
%   A = [_,   theta_0, theta_1, theta_2, ... theta_N,
%        V_0, Y_00,    Y_01,    Y_02,    ... Y_0N,
%        V_1, Y_01,    Y_11,    Y_12,    ... Y_1N,
%        V_2, Y_02,    Y_21,    Y_22,    ... Y_2N,
%         :    :        :        :            :
%        V_M, Y_M0,    Y_M1,    Y_M2,    ... Y_MN    ]
%
%   theta_0 to theta_N are the N+1 evenly-spaced steering angles. theta_0 must
%   be equal to 0. Internally, the VCU mirrors this axis when looking up
%   negative steering angles. Steering angles greater than theta_N are clamped
%   to theta_N. These values are in degrees.
%
%   V_0 to V_M are the M+1 evenly-spaced speed values. V_0 must be equal to 0.
%   These values are in meters per second (converted internally to km/h).
%
%   Y_11 to Y_MN are the M*N yaw-rate values. The Y_0x and Y_x0 values must all
%   be 0. These values are in degrees per second.
%
%   _ is ignored.
%
function x = vcu_Yaw_table(A)

  B = A(2:end, 2:end);

  width = size (B, 2);
  height = size (B, 1);

  if (A (1, 2) > 0.0001 || A (1, 2) < -0.0001)
    printf ("theta_0 must equal 0. Use 'help vcu_yaw_table' for more details.\r\n")
    return
  end

  if (A (2, 1) > 0.0001 || A (2, 1) < -0.0001)
    printf ("V_0 must equal 0. Use 'help vcu_yaw_table' for more details.\r\n")
    return
  end

  for (x = 1:width)
    if (B (1, x) > 0.0001 || B (1, x) < -0.0001)
      printf ("Y_0x values must all equal 0. Use 'help vcu_yaw_table' for more details.\r\n")
      return
    end
  end

  for (y = 1:height)
    if (B (y, 1) > 0.0001 || B (y, 1) < -0.0001)
      printf ("Y_x0 values must all equal 0. Use 'help vcu_yaw_table' for more details.\r\n")
      return
    end
  end

  printf ("\r\n-- BEGIN C SNIPPET --\r\n\r\n");

  printf ("/// @brief The maximum steering angle defined by the lookup table, in degrees.\r\n");
  printf ("#define STEERING_ANGLE_MAX %.4f\r\n", A (1, end));
  printf ("\r\n");
  printf ("/// @brief The width of the lookup table along the steering angle axis. Note this axis is mirrored.\r\n");
  printf ("#define STEERING_ANGLE_WIDTH %i\r\n", width);
  printf ("\r\n")
  printf ("/// @brief The maximum vehicle speed defined by the lookup table, in km/h.\r\n");
  printf ("#define VEHICLE_SPEED_MAX %.4f\r\n", A (end, 1) * 3.6);
  printf ("\r\n")
  printf ("/// @brief The width of the lookup table along the vehicle speed axis. This axis is not mirrored.\r\n");
  printf ("#define VEHICLE_SPEED_WIDTH %i\r\n", height);
  printf ("\r\n");
  printf ("/// @brief The lookup table of ideal yaw rates, in degrees per second.\r\n");
  printf ("static const float YAW_LOOKUP_TABLE [VEHICLE_SPEED_WIDTH][STEERING_ANGLE_WIDTH] =\r\n");
  printf ("{\r\n");

  for (y = 1:height)

    printf ("\t{ ");

    for (x = 1:width)

      printf ("%8.4f", B (y, x));

      if (x ~= width)
        printf (", ");
      else
        printf (" }");
      end

    end

    if (y ~= height)
      printf (",\r\n");
    else
      printf ("\r\n");
    end

  end

  printf ("};\r\n");

  printf ("\r\n-- END C SNIPPET --\r\n\r\n")
    
end

printYawTable (A)
