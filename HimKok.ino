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
  B00100,
  B00100,
  B00100,
  B01110,
  B00100,
  B00100,
  B00100,
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
bool textHelp;
int itextHelp;

// Step mashing parameters and variables
bool stepOn;
int totalStep = 4;
int currentStep = 1;
long stepArr[10][2];
String tempValue;


int inInt;
String mode = "stop";
bool logMode = false;


void setup() {
  lcd.begin(16, 2);                             // Start up of LCD display
  lcd.createChar(0, upArrow);
  lcd.createChar(1, downArrow);
  lcd.createChar(2, line);
  writeLCD("Ikke faa panikk!", "Vennligst vent!", false);
  delay(1000);

  pinMode(heater, OUTPUT);                      // Assigning relay output pin
  digitalWrite(heater, LOW);                    // Turn off Relay (Ensuring it is of
  Serial.begin(9600);                           // Start Serial for debug options
  
  startTemperaturSensor();                      // Start temperatur sensor

  stepMashingSetup();
  
  writeLCD("Brew cntrl ready", "Happy brewin'!!!", false);
  delay(500);
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
    itextHelp ++;
    if (itextHelp == 2) {
      textHelp = !textHelp;
      itextHelp = 0;
    }
    
    // Pre-treatment
    mv = max(min (mv, 100.0), 0.0);
    sp = max(min (sp, 100.0), 0.0);
    err = sp - mv;
    dt_live = millis() - time_1;
    time_1 = millis();
    unsigned long delayAdjust = max(0, dt - dt_live);

    // Timer logic
    if (timerOn && abs(mv-sp)<hyst){
      if (countDown){
        timerTotal -= dt_live;
      } else {
      timerTotal += dt_live;
      }
    }
    if (abs(mv-sp)<hyst && mode == "auto"){
      timerOn = true;
    } else {
      timerOn = false;
    }
    if (stepOn && timerTotal < 0 && currentStep < totalStep-1){
      currentStep ++;
      sp = stepArr[currentStep][0];
      timerTotal = stepArr[currentStep][1];
    } else if (stepOn && timerTotal < 0 ){
      mode = "stop";
      currentStep ++;
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
    if (stepOn && timerTotal < 0) stepOn = false;
  }
  
  if ((lastKeyPadStroke < (millis() - lingerKeyPad)) && inString != "") { // clear inString memory after lingerTime
    inString = "";
  }
  
  if (inString.length() > 0 && customKey ) {
    inString +=  customKey;
  }
  
  if (customKey == '*') {
    inString = customKey;
  }
  if (inString.startsWith("*9")){
    inString = "Out:";
  }
  if (inString.startsWith("*8")){
    inString = "SP:";
  }
  if (inString.startsWith("*7")){
    inString = "Timer:";
  }
  
  if (inString.endsWith("#")) {
    if (inString.startsWith("*")) inString.remove(0, 1);
    inString.remove(inString.length() - 1, 1);
    inInt = inString.toInt();

    if (inString.startsWith("Out:") && inString.length() > 4 ) {
      inString.remove(0, 4);
      out = inString.toFloat();
    } else if (inString.startsWith("SP:") && inString.length() > 3 ) {
      inString.remove(0, 3);
      sp = inString.toFloat();
    } else if (inString.startsWith("Timer:") && inString.length() > 6 ) {
      inString.remove(0, 6);
      timerTotal = inString.toInt()*60000;
      countDown = true;
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
    inString = "";
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
      lcd.write(byte(1));
    } else if (timerOn) {
      lcd.write(byte(0));
    } else {
      lcd.write(byte(2));
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
      if (textHelp){
        firstLine = "0stop 1Man 2Auto";
      } else {
        firstLine = "7Time 8SP  9Out ";
      }
   
      secondLine = inString;
      if (secondLine.startsWith("*") && secondLine.length() > 1){
        secondLine.remove(0,1);
      }
    } else {
       
      timerMin = (timerTotal/60000);
      if (timerMin == 0){
        timerSec = ((timerTotal-(timerMin*60000))/1000);
      } else {
        timerSec = abs((timerTotal-(timerMin*60000))/1000);
      }
      if (stepOn && timerTotal < 0){
        secondLine = "Done!";
      } else{
        secondLine = String(timerMin) + "m" + String(timerSec) + "s";
      }
    
    }
    while ((secondLine.length() < 7 && inString.length() == 0) || (secondLine.length() < 8 && inString.length() != 0) ) {
      secondLine += " ";
    }
    secondLine += mode;
    while ((secondLine.length() < 12 && inString.length() == 0) || (secondLine.length() < 13 && inString.length() != 0) ){
      secondLine += " ";
    }
    secondLine += "H:";
    secondLine += digitalRead(heater);
    
    writeLCD(firstLine, secondLine, (inString.length() == 0));

    lcdNextUpdate = millis() + lcdDelayUpdate;
  }
}

void stepMashingSetup(){
  writeLCD("Steg mesking?   ", "1=Ja og 0=nei   ", false);
  while (inString == "") {
    customKey = customKeypad.getKey();
    inString = customKey;
    if (inString == "1") {                                    // record when last keypress was performed
      stepOn = true;
      Serial.println(customKey);
    } else if (inString == "0"){
      stepOn = false;
    } else {
      inString = "";
    }
  }
  inString = "";

  if (stepOn){
    while (inString != "#") {
      writeLCD("Hvor mange steg?", "Antall steg:" + String(totalStep) + "   ", false);
      customKey = customKeypad.getKey();
      inString = customKey;
      if (inString != "*" && inString != "#" && inString != "0" && inString.length() > 0) totalStep = inString.toInt(); 
    }
    inString = "";
    int i;
    for (i = 0; i < totalStep; i++) {
      // Temperatur setup for the step
      while (inString != "#") {
        writeLCD("Steg" + String(i+1) + " temperatur", "Temperatur:" + tempValue + "C    ", false);
        customKey = customKeypad.getKey();
        inString = customKey;
        if (inString == "*" ){
          tempValue = "";
        } else if ( inString != "#"  && inString.length() > 0) {
          tempValue += inString;
          if (tempValue.length() > 2) tempValue.remove(0,1);
        }
      }
      Serial.println(i);
      stepArr[i][0] = tempValue.toInt();
      tempValue = "";
      inString = "";
      // Timer setup for the step
      while (inString != "#") {
        writeLCD("Steg" + String(i+1) + " tidslengde", "Tid:" + tempValue + "min         ", false);
        customKey = customKeypad.getKey();
        inString = customKey;
        if (inString == "*" ){
          tempValue = "";
        } else if ( inString != "#"  && inString.length() > 0) {
          tempValue += inString;
          if (tempValue.length() > 3) tempValue.remove(0,1);
        }
      }
      stepArr[i][1] = tempValue.toInt()*60000;
      tempValue = "";
      inString = "";
    }
    while (inString != "#") {
      writeLCD("Steg mesk klar! ", "Trykk# for start", false);
      customKey = customKeypad.getKey();
      inString = customKey;
    }
    countDown = true;
    mode = "auto";
    sp = stepArr[0][0];
    timerTotal = stepArr[0][1];
    currentStep = 0;
  }
}
