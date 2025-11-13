# PLEASE DO NOT DISTRIBUTE
# PROPERTY OF EMILY PHILLIPS AND JESSICA SCHLARBAUM

# Import required Python libraries
import RPi.GPIO as GPIO
import time
from gpiozero import MotionSensor
from datetime import datetime
import picamera
import os


#MODIFY IF NEEDED
#23 is the connection port the motion sensor is connected to. Threshold is amount of "motion" for it to be triggered
pir = MotionSensor(23, threshold=.99, queue_len=25)
print(pir)


# MODIFY IF NEEDED. Define GPIO to use on Pi
# This is the connection port that you had connected to the motion sensor to
GPIO_PIR = 23

## MODIFY IF NEEDED. Set the hours that you want. 600 refers to 6:00am and 23:40 refers to 11:40pm
hour_start_1 = 500
hour_end_1 = 2130

print("starting camera")
print(int(str(datetime.now().hour) + str(datetime.now().minute).zfill(2)))  # print current time
# Use BCM GPIO references
# instead of physical pin numbers
GPIO.setmode(GPIO.BCM)

print("PIR Module Test (CTRL-C to exit)")
# Set pin as input and set states
GPIO.setup(GPIO_PIR, GPIO.IN)
Current_State = 0
Previous_State = 0


# Loop until PIR output is 0
print("Waiting for PIR to settle ..")
while GPIO.input(GPIO_PIR) == 1:
  Current_State = 0

# Inifitine loop to run once set up. This loop is split into two pieces, which transfer depending on
# if the current time is within the range secified above
while True:
  try:
    # Check that the current time is within the range
    while int(str(datetime.now().hour) + str(datetime.now().minute).zfill(2)) > hour_start_1 and int(str(datetime.now().hour) + str(datetime.now().minute).zfill(2)) < hour_end_1:

      # Turn on the camera
      with picamera.PiCamera() as camera:
        try:
          time.sleep(3)
          while True:
            # MODIFY IF NEEDED. Change the camera resolution and framerate to get higher quality, larger size videos. Google PiCam for more setting options
            camera.resolution = (640, 480)
            camera.framerate = 60
            camera.annotate_background = picamera.Color('black')
            # stream seconds is the buffer seconds (before motion detected)
            stream = picamera.PiCameraCircularIO(camera, seconds=5)
            # MODIFY IF NEEDED. You can adjust the camera recording and this will change how large the saved files are
            camera.start_recording(stream, format='h264', quality=25)

            print("streaming")

            # While loop to check for motion
            while True:
              print("looking for motion")
              pir.wait_for_motion(10)

              if pir.motion_detected:
                print("motion detected!")
                # camera.wait seconds indiates how much video to take AFTER buffer
                camera.annotate_text = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                camera.wait_recording(10)
                timestamp = datetime.now().strftime('%d-%m-%y_%H-%M-%S')

                # MODIFY IF NEEDED. Create the file path for where to save the video. Your file path depends on your raspberry pi.
                outputvideo = '/mnt/usb/videos/{}.h264'.format(timestamp)
                # Mount the USB, save it to the usb, then unmount it. This prevents the files from corrupting
                # when the battery source dies with the usb mounted
                print("mounting usb")
                os.system("sudo mkdir /mnt/usb/videos")
                os.system("sudo mount /dev/sda1 /mnt/usb/videos")

                # stream.copy_to seconds indicates how many seconds of video to save (buffer seconds + post-buffer seconds). without i$
                stream.copy_to(outputvideo, seconds=15)
                print("unmounting")
                os.system("sudo umount /mnt/usb/videos")
                print(outputvideo)
                print("copied to usb")
                print("sleeping")
                time.sleep(3)

              # Check that the time is still valid. If not, break the loop
              if int(str(datetime.now().hour) + str(datetime.now().minute).zfill(2)) < hour_start_1 or int(str(datetime.now().hour) + str(datetime.now().minute).zfill(2)) > hour_end_1:
                print("time for sleep now")
                break

            break

        # Allows you to quit the program if you need to. Use CTRL-C
        except KeyboardInterrupt:
          print("Quit")
          camera.stop_recording()

    # Printing sleep
    print("zzzzzzzz")

# Allows you to quit the program if you need to. Use CTRL-C
  except KeyboardInterrupt:
   print("Quit")
   time.sleep(5)



