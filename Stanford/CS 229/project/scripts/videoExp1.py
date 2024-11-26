import argparse
import os
from ai.dive.data.label_reader import LabelReader
from ai.dive.data.image_file_classification import ImageFileClassificationDataset
from transformers import ViTForImageClassification, ViTImageProcessor
from transformers import TrainingArguments
from transformers import Trainer
from datasets import load_metric, load_dataset
import torch
import numpy as np

# We must take the dataset and stack it into pytorch tensors
# Our batch size above was 16
# So this will be a stack of 16 images into a tensor
def collate_fn(batch):
    return {
        'pixel_values': torch.stack([x['pixel_values'] for x in batch]),
        'labels': torch.tensor([x['labels'] for x in batch])
    }

# We want to evaluate accuracy of the model on the test/eval set
metric = load_metric("accuracy")
def compute_metrics(p):
    return metric.compute(predictions=np.argmax(p.predictions, axis=1), references=p.label_ids)

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Train a ViT on dataset')
    parser.add_argument('-d', '--data', default="../datasets/DAiSEE/Labels", type=str, help='datasets to train/eval model on')
    parser.add_argument('-o', '--output', default="../trained_models", type=str, help='output file to write results to')
    parser.add_argument('-m', '--base_model', default="google/vit-base-patch16-224-in21k", type=str, help='The base model to use')
    parser.add_argument('-g', '--gpu', default=True, help='Train on the GPU if supported')
    args = parser.parse_args()

    # Parse Args
    # ...
    labels_file = os.path.join(args.data, "labels.txt")
    label_reader = LabelReader(labels_file)
    labels = label_reader.labels()
    print(labels)
    # ....

    # Same processor as before that is tied to the model
    processor = ViTImageProcessor.from_pretrained(args.base_model)

    # Load the dataset into memory, and convert to a hugging face dataset
    print("Preparing train dataset...")
    train_file = "test.csv"    ## modified
    ds = ImageFileClassificationDataset(
        data_dir=args.data,
        file=train_file,
            label_reader=label_reader,
        img_processor=processor,
    )
    train_dataset = ds.to_hf_dataset()
    
    test_file = "validation.csv"
    ds = ImageFileClassificationDataset(
        data_dir=args.data,
        file=test_file,
        label_reader=label_reader,
        img_processor=processor
    )
    eval_dataset = ds.to_hf_dataset()

    print(train_dataset[0])
    print(train_dataset[0]['pixel_values'].shape)

    model = ViTForImageClassification.from_pretrained(
        args.base_model,
        num_labels=len(labels),
        id2label={str(i): c for i, c in enumerate(labels)},
        label2id={c: str(i) for i, c in enumerate(labels)}
    )

    training_args = TrainingArguments(
        output_dir=args.output, # directory to save the model
        per_device_train_batch_size=16,
        evaluation_strategy="steps",
        num_train_epochs=4, # loop through the data N times
        no_cuda=(not args.gpu), # use the GPU or not
        save_steps=1000, # save the model every N steps
        eval_steps=1000, # evaluate the model every N steps
        logging_steps=10,
        learning_rate=2e-4,
        save_total_limit=2, # only keep the last N models
        remove_unused_columns=False,
        report_to='tensorboard',
        load_best_model_at_end=True,
    )
    # Instantiate a trainer with all the components we have built so far
    trainer = Trainer(
        model=model,
        args=training_args,
        data_collator=collate_fn,
        compute_metrics=compute_metrics,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset,
        tokenizer=processor,
    )

    # Kick off the train
    print("Training model...")
    train_results = trainer.train()
    trainer.save_model()
    trainer.log_metrics("train", train_results.metrics)
    trainer.save_metrics("train", train_results.metrics)
    trainer.save_state()

if __name__ == '__main__':
    main()