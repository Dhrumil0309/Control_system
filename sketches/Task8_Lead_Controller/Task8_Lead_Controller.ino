// ============================================================================
//  Task 8: Lead Compensator – Hardware Implementation & Tuning
//  ES 245 Course Project Part-II  |  Ball/Car-and-Beam System
// ============================================================================
//
//  Plant (Task 5):   G(s) = -0.6463 / [ s (s + 0.0170) ]   [cm / deg]
//  Compensator (Task 7):
//      C(s) = 4.50 * (s + 0.050) / (s + 3.50)       [deg / cm]
//
//  Discrete C(z) via Tustin at Ts = 0.02 s (50 Hz):
//      C(z) = (4.350000 z - 4.345652) / (z - 0.932367)
//
//  Difference equation:
//      u[k] = 0.932367*u[k-1] + 4.350000*e[k] - 4.345652*e[k-1]
//      servo_deg = u0_trim + DIRECTION * u[k]
//
//  Sensor: HC-SR04 Ultrasonic (trigPin=9, echoPin=A1)
//  Actuator: Standard hobby servo on pin 3
//
//  Serial Tuning Commands (send via Serial Monitor at 115200 baud):
//      G<val>   Set gain multiplier   (e.g. G1.2 → scale controller by 1.2x)
//      T<val>   Set servo trim angle  (e.g. T108.0)
//      R<val>   Set setpoint distance (e.g. R25.0 → target 25 cm)
//      D        Flip servo direction  (toggles DIRECTION between +1 / -1)
//      S        Start / Stop control  (toggle on/off)
//      ?        Print current settings
//
// ============================================================================

#include <Servo.h>

// ─────────────────────────── Pin Definitions ────────────────────────────────
const int TRIG_PIN  = 9;
const int ECHO_PIN  = A1;
const int SERVO_PIN = 3;

// ─────────────────────────── Hardware Configuration ─────────────────────────
float u0_trim   = 107.6;     // Servo neutral trim (deg) – identified in Task 5
float setpoint  = 20.0;      // Target distance (cm)

const float SERVO_MIN = 50.0;   // Mechanical lower safety limit (deg)
const float SERVO_MAX = 160.0;  // Mechanical upper safety limit (deg)
int   direction = -1;           // -1 because plant gain K is negative
                                // (increasing servo angle → decreasing distance)
                                // Flip to +1 if your physical setup is inverted

// ─────────────────── Discrete Lead Compensator Coefficients ─────────────────
// C(z) = (b0*z + b1) / (z + a1)
//      = (4.350000*z - 4.345652) / (z - 0.932367)
//
// From C(s) = 4.50 * (s + 0.050) / (s + 3.50) via bilinear transform at Ts=0.02s
//
// Difference eq:  u[k] = (-a1)*u[k-1] + b0*e[k] + b1*e[k-1]

const float B0_BASE =  4.350000;   // Numerator z^1 coefficient
const float B1_BASE = -4.345652;   // Numerator z^0 coefficient
const float A1_NEG  =  0.932367;   // = -a1 (denominator z^0, sign-flipped)

float gainMultiplier = 1.0;        // Live gain scaling factor for tuning

// ─────────────────────── Controller State Variables ─────────────────────────
float u_prev     = 0.0;    // Previous controller output u[k-1]
float error_prev = 0.0;    // Previous error e[k-1]

// ─────────────────────── Sensor Filtering ───────────────────────────────────
const float ALPHA_DIST   = 0.30;   // EMA filter weight (higher = less smoothing)
const float MAX_JUMP_CM  = 8.0;    // Reject readings that jump more than this in one step
float filteredDist       = 20.0;   // Filtered distance (cm)
float lastValidDist      = 20.0;   // Last accepted raw reading (cm)

// ─────────────────────── Timing ─────────────────────────────────────────────
const unsigned long SAMPLE_TIME_MS = 20;   // 50 Hz control loop (Ts = 0.02 s)
unsigned long lastSampleTime = 0;

// ─────────────────────── Control State ──────────────────────────────────────
bool controlEnabled = false;   // Starts DISABLED – enable with 'S' command
                               // or wait for auto-start after 3 s

Servo servo;

// ============================================================================
//  Ultrasonic Sensor: Single-shot distance measurement (cm)
//  Returns -1.0 on timeout / out-of-range
// ============================================================================
float readUltrasonicCm() {
  // Trigger pulse
  digitalWrite(TRIG_PIN, LOW);
  delayMicroseconds(2);
  digitalWrite(TRIG_PIN, HIGH);
  delayMicroseconds(10);
  digitalWrite(TRIG_PIN, LOW);

  // Measure echo duration (timeout 25 ms ≈ 430 cm max)
  long duration = pulseIn(ECHO_PIN, HIGH, 25000UL);

  if (duration == 0) return -1.0;   // No echo received

  float dist_cm = duration * 0.0343 / 2.0;

  // Reject readings outside useful sensor range
  if (dist_cm < 2.0 || dist_cm > 50.0) return -1.0;

  return dist_cm;
}

// ============================================================================
//  Filtered distance reading with outlier rejection + EMA
// ============================================================================
float getFilteredDistance() {
  float raw = readUltrasonicCm();

  // If sensor returned invalid, hold previous value
  if (raw < 0.0) return filteredDist;

  // Reject sudden jumps (noise spikes)
  if (abs(raw - lastValidDist) > MAX_JUMP_CM) return filteredDist;

  lastValidDist = raw;

  // Exponential moving average
  filteredDist = ALPHA_DIST * raw + (1.0 - ALPHA_DIST) * filteredDist;
  return filteredDist;
}

// ============================================================================
//  Parse serial tuning commands
// ============================================================================
void processSerialCommand() {
  if (!Serial.available()) return;

  String cmd = Serial.readStringUntil('\n');
  cmd.trim();
  if (cmd.length() == 0) return;

  char tag = cmd.charAt(0);
  float val = 0.0;
  if (cmd.length() > 1) {
    val = cmd.substring(1).toFloat();
  }

  switch (tag) {
    case 'G': case 'g':   // Gain multiplier
      if (val > 0.0 && val < 20.0) {
        gainMultiplier = val;
        Serial.print(">> Gain multiplier set to "); Serial.println(gainMultiplier, 4);
      }
      break;

    case 'T': case 't':   // Trim angle
      if (val >= 60.0 && val <= 150.0) {
        u0_trim = val;
        Serial.print(">> Trim angle set to "); Serial.print(u0_trim, 1); Serial.println(" deg");
      }
      break;

    case 'R': case 'r':   // Setpoint distance
      if (val >= 5.0 && val <= 45.0) {
        setpoint = val;
        // Reset controller state for clean transient
        u_prev = 0.0;
        error_prev = 0.0;
        Serial.print(">> Setpoint set to "); Serial.print(setpoint, 1); Serial.println(" cm");
      }
      break;

    case 'D': case 'd':   // Flip direction
      direction = -direction;
      Serial.print(">> Direction flipped to "); Serial.println(direction);
      break;

    case 'S': case 's':   // Start / Stop toggle
      controlEnabled = !controlEnabled;
      if (controlEnabled) {
        // Reset controller state on re-enable
        u_prev = 0.0;
        error_prev = 0.0;
        Serial.println(">> CONTROL ENABLED");
      } else {
        servo.write((int)u0_trim);  // Return to neutral
        Serial.println(">> CONTROL DISABLED (servo at trim)");
      }
      break;

    case '?':   // Print current settings
      Serial.println("─── Current Settings ───");
      Serial.print("  Setpoint:    "); Serial.print(setpoint, 1);     Serial.println(" cm");
      Serial.print("  Trim (u0):   "); Serial.print(u0_trim, 1);     Serial.println(" deg");
      Serial.print("  Gain mult:   "); Serial.println(gainMultiplier, 4);
      Serial.print("  Direction:   "); Serial.println(direction);
      Serial.print("  Control:     "); Serial.println(controlEnabled ? "ENABLED" : "DISABLED");
      Serial.print("  Distance:    "); Serial.print(filteredDist, 1); Serial.println(" cm");
      Serial.println("────────────────────────");
      break;

    default:
      Serial.println(">> Unknown command. Use G/T/R/D/S/?");
      break;
  }
}

// ============================================================================
//  SETUP
// ============================================================================
void setup() {
  Serial.begin(115200);

  pinMode(TRIG_PIN, OUTPUT);
  pinMode(ECHO_PIN, INPUT);

  servo.attach(SERVO_PIN);
  servo.write((int)u0_trim);   // Start at neutral trim

  Serial.println();
  Serial.println("============================================");
  Serial.println("  Task 8: Lead Compensator Controller");
  Serial.println("  C(s) = 4.50*(s+0.05)/(s+3.50)");
  Serial.println("============================================");
  Serial.println("Send '?' for settings, 'S' to start control");
  Serial.println("Serial Plotter columns:");
  Serial.println("  Setpoint_cm, Distance_cm, Servo_deg, Error_cm");
  Serial.println("============================================");

  // Take a few readings to initialise the filter
  delay(500);
  for (int i = 0; i < 10; i++) {
    float d = readUltrasonicCm();
    if (d > 0) {
      filteredDist = d;
      lastValidDist = d;
    }
    delay(25);
  }

  Serial.print("Initial distance: ");
  Serial.print(filteredDist, 1);
  Serial.println(" cm");
  Serial.println("Waiting 3 s at trim, then auto-enabling control...");

  delay(3000);    // Let the beam settle at trim for 3 seconds

  // Auto-enable control after startup delay
  controlEnabled = true;
  u_prev = 0.0;
  error_prev = 0.0;
  lastSampleTime = millis();

  Serial.println(">> CONTROL ENABLED (auto-start)");
}

// ============================================================================
//  MAIN LOOP – 50 Hz Discrete Lead Compensator
// ============================================================================
void loop() {
  // ── Handle serial tuning commands (non-blocking) ──
  processSerialCommand();

  // ── Enforce fixed 50 Hz sampling rate ──
  unsigned long now = millis();
  if (now - lastSampleTime < SAMPLE_TIME_MS) return;
  lastSampleTime = now;

  // ── 1. Read and filter sensor distance ──
  float dist_cm = getFilteredDistance();

  // ── 2. Compute control output ──
  float servo_deg;

  if (controlEnabled) {
    // Error: e[k] = setpoint - measured  (positive = car too close to sensor)
    float error = setpoint - dist_cm;

    // Discrete Lead Compensator difference equation:
    //   u[k] = A1_NEG * u[k-1]  +  gainMult * B0 * e[k]  +  gainMult * B1 * e[k-1]
    //
    // gainMultiplier scales the overall compensator gain for live tuning.
    // Default gainMultiplier = 1.0 gives the exact Task 7 design.
    float u = A1_NEG * u_prev
            + gainMultiplier * B0_BASE * error
            + gainMultiplier * B1_BASE * error_prev;

    // Update state for next iteration
    u_prev     = u;
    error_prev = error;

    // Servo command: trim + direction * controller output
    // direction = -1 compensates for negative plant gain K = -0.6463
    servo_deg = u0_trim + (float)direction * u;

    // Mechanical saturation
    servo_deg = constrain(servo_deg, SERVO_MIN, SERVO_MAX);

    servo.write((int)servo_deg);

  } else {
    // Control disabled: hold at trim
    servo_deg = u0_trim;
    servo.write((int)servo_deg);
  }

  // ── 3. Stream data for Serial Plotter / CSV logging ──
  // Format: Setpoint_cm, Distance_cm, Servo_deg, Error_cm
  Serial.print("Setpoint_cm:");  Serial.print(setpoint, 1);            Serial.print(",");
  Serial.print("Distance_cm:");  Serial.print(dist_cm, 1);             Serial.print(",");
  Serial.print("Servo_deg:");    Serial.print(servo_deg, 1);           Serial.print(",");
  Serial.print("Error_cm:");     Serial.println(setpoint - dist_cm, 2);
}
