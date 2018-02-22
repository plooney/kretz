# README #

### What is this repository for? ###

* OXCNN_core is a barebones python program designed to be extended to create more complex neural networks for patch based image analysis
* Version 0.01

### How do I get set up? ###

**Create a set of test cases** 

`python3 -m oxcnn.test.utils ~/Desktop/TestVolumes`

**Write out the TensorFlowRecords**

`python3 main.py --model models.deepmedic write --save_dir ~/Desktop/TestRect-tfr --data_dir ~/Desktop/TestVolumes/`

**Train the model**

`python3 main.py --model models.deepmedic train --tfr_dir ~/Desktop/TestRect-tfr --save_dir ~/Desktop/TestRect-out --test_data ~/Desktop/TestRect-tfr/meta_data.txt --num_epochs 10 --batch_size 10`

**Test the model**

`python3 main.py --model models.deepmedic test --save_dir ~/Desktop/TestRect-out/test --test_data_file ~/Desktop/TestRect-tfr/meta_data.txt --model_file /home/padraig/Desktop/TestRect-out/epoch_model.ckpt-2243.meta --batch_size 10`

**Dependencies**

`tflearn (v0.3), tensorflow (v1.2), pandas, nibabel`

Python 3.5 recommended.

** How to run tests **

`python3 -m unittest`
