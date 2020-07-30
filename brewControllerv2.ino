/*
   Program: BrewController
   Version: 1.2
   Author: Kristian Dyb Strand
   Date: 22.12.19


*/
#include <Wire.h>
#include <OneWire.h>            // used by temperature sensor
#include <DallasTemperature.h> 
#include <LiquidCrystal.h>
#include <Keypad.h>

// Setting up Keypad
const byte ROWS = 4; //four rows
const byte COLS = 3; //four columns

char hexaKeys[ROWS][COLS] = {
  {'1', '2', '3'},
  {'4', '5', '6'},
  {'7', '8', '9'},
  {'*', '0', '#'}
};                                      //define the symbols on the buttons of the keypads
byte rowPins[ROWS] = {A0, A1, A2, A3};  // Data pins for row
byte colPins[COLS] = {A4, A5, 3};       // Data pins for collumn 
Keypad customKeypad = Keypad( makeKeymap(hexaKeys), rowPins, colPins, ROWS, COLS);  // Assigning Keypad variable
unsigned long lastKeyPadStroke;         // Assigning last input variable
unsigned long lingerKeyPad = 10000;     // Assigning variable for duration of keypad memory
String inString = "";                   // Assigning variable to collect input
char customKey;                         // Collection variable for keystrokes

// Setting up Temperature probe (DS18B20)
OneWire ds(5);                                        // Data pin for sensor
DallasTemperature sensors(&ds);                       // Assigning sensor variable
DeviceAddress sensorAddress;                          // Assinging sensor address variable
int resolution = 12;                                  // Setting sensor resolution (8, 10, 11, 12 bits)
unsigned long lastTempRequest = 0;                    // Initialise timestamp used to controll update rate of sensor
int  delaySensor =  750 / (1 << (12 - resolution)); // Initialise sensor update rate based on resolution

// Setting up LCD screen
const int rs = 12, en = 11, d4 = 10, d5 = 9, d6 = 8, d7 = 7;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);
unsigned long lcdNextUpdate;
unsigned long lcdDelayUpdate = 500;
byte upArrow[8] = {
  B00100,
  B01110,
  B10101,
  B00100,
  B00100,
  B00100,
  B00100,
};
byte downArrow[8] = {
  B00100,
  B00100,
  B00100,
  B00100,
  B10101,
  B01110,
  B00100,
};
byte line[8] = {
  B00000,
  B00000,
  B00000,
  B11111,
  B00000,
  B00000,
  B00000,
};



// Setting up Heater control
const int heater = 13;
float pulseLength = 10000;     // Pulse duration in milliseconds


// PID parameters
float kp = 73.2;        // Propotional band 17.4
float ti = 240000;    // Integral time
float dt = 100;       //total cycle time (ms)

// Variables
float out = 100;      // Output value
float err = 0;        // Error (MV - SP)
float err_1 = 0;      // Last error
float mv = 0;         // measuring value
float sp = 30;        //Setpoint
bool hystOn;
float hyst = 0.5;

// Filtering parameters and variables
const int numReadings = 10;
float readings[numReadings];      // the readings from the analog input
int readIndex = 0;              // the index of the current reading
float total = 0;                  // the running total
//float average = 0;                // the average

unsigned long pulseStart;
unsigned long time_1;
unsigned long dt_live;

long timerTotal;
long timerSec;
long timerMin;
bool timerOn;
bool countDown;


int inInt;
String mode = "stop";
bool logMode = false;


void setup() {
  lcd.begin(16, 2);                             // Start up of LCD display
  writeLCD("Brew cntrl init", "Wait for it!", false);
  delay(2000);
  lcd.createChar(0, upArrow);
  lcd.createChar(1, downArrow);
  lcd.createChar(2, line);

  writeLCD("Init heater pin", "Wait for it!", false);
  pinMode(heater, OUTPUT);                      // Assigning relay output pin
  digitalWrite(heater, LOW);                    // Turn off Relay (Ensuring it is of
  delay(2000);
  Serial.begin(9600);
  
  writeLCD("Start tempsensor", "Wait for it!", false);
  startTemperaturSensor();                      // Start temperatur sensor

  
  delay(2000);
  
  writeLCD("Brew cntrl ready", "Happy brewin'!!!", false);
  delay(2000);
  updateLCD;
  pulseStart = millis();                         // Initialising the pulse start
  lcdNextUpdate = millis() + lcdDelayUpdate;     // Initialising the lcd update cycle
}

void(* resetFunc) (void) = 0;//declare reset function at address 0

void loop() {
  updateLCD();
  keyPadRead();

  if (millis() - lastTempRequest >= delaySensor){
    updateMeasurement();
  
    // Pre-treatment
    mv = max(min (mv, 100.0), 0.0);
    sp = max(min (sp, 100.0), 0.0);
    err = sp - mv;
    dt_live = millis() - time_1;
    time_1 = millis();
    unsigned long delayAdjust = max(0, dt - dt_live);

    // Timer logic
    if (timerOn && abs(mv-sp)<hyst){
      timerTotal =+ dt_live;
      timerMin = (timerTotal/60000);
      timerSec = ((timerTotal-(timerMin*60000))/1000);
    }
    if (abs(mv-sp)<hyst){
      timerOn = true;
    } else {
      timerOn = false;
    }
    
  
    // PI algoritme
    if (mode == "auto") {
      // Proportional band
      out = out + (kp * (err - err_1));
      // Integral band
      out = out + ((kp / ti) * err * float(dt_live));
    } else if (mode == "stop") {
      out = 0.0;
    }  else if (mode == "FP") {
      out = 100.0;
    }
    err_1 = err;
  }
  // Output treatment
  out = constrain(out, 0.0, 100.0);
  
  
  // Hysteris control algoritme
  if ((mv - sp) > 0) {
    hystOn = false;
  } else if ((sp - mv) > hyst) {
    hystOn = true;
  }

  // Pulse modulation
  unsigned long now = millis();
  if ((mode == "hyst" && hystOn) || mode != "hyst") {
    int pulseHeater = (out / 100) * pulseLength;
    if (now - pulseStart > pulseLength) {
      pulseStart = now;
    }
    if (pulseHeater > now - pulseStart) {
      digitalWrite(heater, HIGH);
    } else {
      digitalWrite(heater, LOW);
    }
  } else {
    digitalWrite(heater, LOW);
  }

}



void keyPadRead() {

// Function to register keypresses and decoding it to actions
  
  customKey = customKeypad.getKey();
  if (customKey) {                                    // record when last keypress was performed
    lastKeyPadStroke = millis();
  }
  if ((lastKeyPadStroke < (millis() - lingerKeyPad)) && inString != "") { // clear inString memory after lingerTime
    inString = "";

    Serial.println(lingerKeyPad);
  }
  if (inString.startsWith("*") && customKey ) {
    inString +=  customKey;
  }
  if (customKey == '*') {
    inString = customKey;
  }
  if (inString.startsWith("*0")){
    inString = "*Out:";
  }
  if (inString.startsWith("*9")){
    inString = "*SP:";
  }
  
  if (inString.endsWith("#")) {
    inString.remove(0, 1);
    inString.remove(inString.length() - 1, 1);
    inInt = inString.toInt();

    if (inString.startsWith("Out:") && inString.length() > 4 ) {
      inString.remove(0, 4);
      out = inString.toFloat();
    } else if (inString.startsWith("SP:") && inString.length() > 3 ) {
      inString.remove(0, 3);
      sp = inString.toFloat();
    } else {
      switch (inInt) {
        case 0:
          mode = "stop";
          break;
        case 1:
          mode = "man";
          break;
        case 2:
          mode = "auto";
          break;
        case 3:
          mode = "hyst";
          break;
        case 58008:
          mode = "FP";
          break;
        case 1337:
          logMode =  !logMode;
          break;
        case 666:
          resetFunc();
      }
    }
  }
  if (customKey && logMode) {
    Serial.println(inString);
    Serial.println(customKey);
  }

}

void startTemperaturSensor(){
  /*
   Function to start DS18B20 temperatur sensor
   */
  sensors.begin();                                        // Start sensor communication
  sensors.getAddress(sensorAddress, 0);                   // Get sensor address
  sensors.setResolution(sensorAddress, resolution);       // Set sensor resolution (8, 10, 11, 12 bits)
  sensors.requestTemperatures();                          // Poll the sensor to get values first measurment
  mv = sensors.getTempCByIndex(0);                        // Set inital temperatur measurment
  sp = mv;                                                // Initialise setpoint equal to temperature
  sensors.setWaitForConversion(false);                    // When set to False the request does not wait for reply
  sensors.requestTemperatures();                          // Poll the sensor to start the temperature requests
  lastTempRequest = millis();                             // timestamp to keep track of time since request
  
  for (int i = 0; i < numReadings; i++) {
    readings[i] = mv;                           // Initialising measurment smoothing
    total += mv;
  }
  
}

void updateMeasurement(){
  float ms = sensors.getTempCByIndex(0);
  sensors.requestTemperatures();
  lastTempRequest = millis();
  total = total - readings[readIndex];
  readings[readIndex] = ms;
  total = total + readings[readIndex];
  readIndex = readIndex + 1;
  if (readIndex >= numReadings) {
    readIndex = 0;
  }
  mv = total / float(numReadings);
}

void writeLCD(String firstLine, String secondLine, boolean timerShow){
  /*
  Function to write 2 lines to a 16x2 LCD display
   */
  lcd.setCursor(0, 0);
  lcd.print(firstLine);
  lcd.setCursor(0, 1);
  if (timerShow){
    if (countDown && timerOn){
      lcd.print(byte(1));
    } else if (timerOn) {
      lcd.print(byte(0));
    } else {
      lcd.print(byte(2));
    }
  }
  lcd.print(secondLine);
}

void updateLCD() {
  if (lcdNextUpdate < millis()) {
    if (logMode) {
      String logText =  String(mv) + "," + String(sp) + "," + String(out) + "," + String(dt_live);
      Serial.println(logText);
    }
    // Drive the LCD
    String firstLine;
    String secondLine;
    
    firstLine = (String(mv,1) + "|" + String(sp,1)  + "|" + String(out,1));
    while (firstLine.length() < 16){
      firstLine += " ";
    }
    
    if (inString != "") {
      secondLine = inString;
      if (secondLine.startsWith("*") && secondLine.length() > 1){
        secondLine.remove(0,1);
      }
    } else {
      if (timerTotal = 0) {
        secondLine = "mv!=SP";
      } else {
        secondLine = String(timerMin) + "m" + String(timerSec) + "s";
      }
    }
    while (secondLine.length() < 8){
      secondLine += " ";
    }
    secondLine += mode;
    while (secondLine.length() < 13){
      secondLine += " ";
    }
    secondLine += "H:";
    secondLine += digitalRead(heater);

    writeLCD(firstLine, secondLine, true);

    lcdNextUpdate = millis() + lcdDelayUpdate;
  }
}
