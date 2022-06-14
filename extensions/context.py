import os
import uuid
import csv
import re
import unicodedata
import itertools
from copier_templates_extensions import ContextHook


def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


class ContextUpdater(ContextHook):
    def __init__(self, environment, **kwargs):  # noqa: N805 (self)
        """Initialize the object.
        Arguments:
            environment: The Jinja environment.
        """
        super(ContextUpdater, self).__init__(environment, **kwargs)
        
        print(self)
        print(dir(self))
        print(environment)
        print(dir(environment))
        print(dir(ContextUpdater))
        print(kwargs)
        # Generate a unique id that is persisted for all file contexts
        self.unique_id = uuid.uuid4()
        self.working_dir = os.getcwd()

    def get_sample_sheet_path(self, context):
        sample_sheet = context.get("sample_sheet", None)
        template_path = context["_src_path"]
        if not sample_sheet:
            raise FileNotFoundError("You must to specify a sample sheet")
        try:
            # if it's an absolute path already, just try to resolve it and use it
            if os.path.isabs(sample_sheet):
                sample_sheet_path = os.path.realpath(sample_sheet, strict=True)
            else:
                # otherwise, let's try and turn it into an absolute path
                # print(os.path.relpath(sample_sheet, start=template_path))
                sample_sheet_path = os.path.realpath(os.path.join(template_path, sample_sheet), strict=True)
        except:
            raise FileNotFoundError(f"Unable to find sample sheet {sample_sheet}, please make sure you tryped the filename and path properly.")
        return sample_sheet_path

    def read_sample_sheet(self, sample_sheet_path, filter_user_id, filter_experiment_id):
        """Reads the sample sheet file, filters only the samples that match the
        user_id and experiment_id
        """
        samples = []
        try:
            with open(sample_sheet_path) as f:
                reader = csv.DictReader(f, delimiter="\t", fieldnames=["Sample number","User","Sample name","Barcode2","Lane","Barcode","Expt grouping","Reads"], restkey="Reads2")
                headers = next(reader)
                for sample in reader:
                    sample["Reads"] = [read_file.strip() for read_file in [sample["Reads"]] + sample["Reads2"]]
                    expt_grouping = sample.get("Expt grouping", "").strip().split("_")
                    if len(expt_grouping) < 3:
                        raise ValueError(f"Error in sample sheet on line number {reader.line_num}. The experimental grouping specified is invalid. It should be like <user id>_<experiment_id><sample id>_<replicate num>")
                    try:
                        user_id = int(expt_grouping[0])
                        expt_reptype = expt_grouping[1]
                        experiment_id = expt_reptype[0]
                        sample_condition_num = int(expt_reptype[1])
                        replicate_num = int(expt_grouping[2])
                    except:
                        raise ValueError(f"Error in sample sheet on line number {reader.line_num}. The experimental grouping specified is invalid. It should be like <user_id>_<experiment_id><sample_condition_num>_<sample_condition_num>. Probably you are using non-integer values for user_id, sample_condition_num or replicate_num")
                    
                    # skip processing this sample if it doesn't match the specified user id or experiment id
                    if user_id != filter_user_id or experiment_id != filter_experiment_id:
                        next

                    try:
                        sample["Sample number"] = int(sample["Sample number"].strip())
                        sample["User"] = sample["User"].strip()
                        sample["Sample name"] = sample["Sample name"].strip()
                        sample["sample_name_slug"] = slugify(sample["Sample name"])
                        sample["Barcode2"] = sample["Barcode2"].strip()
                    except:
                        raise ValueError(f"There is an error in one of the following columns on line number {reader.line_num}: Sample number, Sample name")

                    # try to pull out the sample lanes and their corresponding read files
                    try:
                        sample["Lane"] = sample["Lane"].strip()
                        if "," in sample["Lane"]:
                            lanes = ["L0{}".format(l) for l in sample["Lane"].split(",")]
                        else:
                            lanes = ["L0{}".format(sample["Lane"])]
                        sample["Lanes"] = lanes
                        # check that there are two read-pairs for every lane specified
                        if len(lanes) * 2 != len(sample["Reads"]):
                            raise ValueError(f"Error in sample sheet on line number {reader.line_num}: you do not have the correct number of read files, there should be two read files per lane specified in the Lanes column.")
                        # grab the read files for each lane specified
                        reads_per_lane = []
                        for i, lane in enumerate(lanes):
                            fq1 = sample["Reads"][i*2]
                            fq2 = sample["Reads"][(i*2)+1]
                            
                            reads_per_lane.append([fq1,fq2])
                        sample["reads_per_lane"] = reads_per_lane

                        # try to create individual units from the samples and reads
                        sample["units"] = []
                        for lane, reads_in_lane in zip(lanes, sample["reads_per_lane"]):
                            sample["units"].append(dict(
                                lane=lane,
                                fq1=reads_in_lane[0],
                                fq2=reads_in_lane[1],
                            ))
                    except:
                        raise ValueError(f"Error in sample sheet on line number {reader.line_num}: You have issues with your lanes or the read file names you specified.")
                    
                    sample["condition_label"] = f"Treatment_{sample_condition_num}" if sample_condition_num > 0 else "Control"
                    sample["user_id"] = user_id
                    sample["experiment_id"] = experiment_id
                    sample["sample_condition_num"] = sample_condition_num
                    sample["replicate_num"] = replicate_num
                    samples.append(sample)
        except:
            raise ValueError(f"There was an error processing the sample sheet {sample_sheet_path}. Please check it for proper formatting. Scroll up to check for additional error messages which may give more detail. Als make sure that all spaces are literal tab characters.")
        
        # prepare samples for template rendering
        # sort them by replicate number, then by treatment condition
        return sorted(sorted([dict(
                sample_num=s["Sample number"],
                id=s["sample_name_slug"],
                barcode=s["Barcode2"],
                condition_label=s["condition_label"],
                # coded condition number
                #   0 = control
                #   1 = experimental condition 1
                #   2 = exp. condition 2, ...
                condition_num=s["sample_condition_num"],
                units=sorted(s["units"], key=lambda u: u["lane"]),
                user_id=s["user_id"],
                experiment_id=s["experiment_id"],
                replicate_num=s["replicate_num"],
            ) for s in samples], key=lambda s: s["replicate_num"]), key=lambda s: s["condition_num"])
        

    def hook(self, context):
        species = context.get("species")
        is_hkust = context.get("is_hkust", False)

        extra_context = {
            "species_friendly": species.title().replace("_", " "),
            "unique_id": self.unique_id,
            "docker_image": f"rnaseq:{self.unique_id}",
        }
        
        
        # template_path = context["_src_path"]
        # # dest_path = context["_copier_conf"]["dst_path"]
        # sample_sheet_path = self.get_sample_sheet_path(context)
        # samples = []
        # comparison_items = []
        # print(self.working_dir)
        # print(dir(self))
        # print(dir(context["_copier_conf"].keys()))
        print(dir(context))
        if is_hkust:
            user_id = context["user_id"]

            # Ensure they specified an experiment id during copier step
            experiment_id = context.get("experiment_id", None)
            if not experiment_id:
                raise ValueError(f"You need to specify an experiment ID indiciating which samples to process.")

            sample_sheet_path = self.get_sample_sheet_path(context)
            extra_context["sample_sheet_path"] = sample_sheet_path

            
            samples = self.read_sample_sheet(sample_sheet_path, user_id, experiment_id)
            extra_context["samples"] = samples

            condition_labels = []
            for condition_label, condition in itertools.groupby(samples, lambda s: s["condition_label"]):
                condition_labels.append(condition_label)
            comparisons = [tuple(reversed(c)) for c in itertools.combinations(condition_labels, 2)]
            comparison_items = [dict(label=f"{c[0]}-vs-{c[1]}", groups=c) for c in comparisons]

            extra_context["diffexp_contrasts"] = comparison_items
            
        return extra_context