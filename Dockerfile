FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.


# Install Roary and R
RUN apt-get update && apt-get install -y roary \
	r-base \
	&& rm -rf /var/lib/apt/lists/*

# Install ggplot2 for graphing functionality
RUN Rscript -e "install.packages('ggplot2', repos = 'http://cran.us.r-project.org')"

# run pip installations
RUN pip install -U pip \
	&& pip install pandas \
	&& pip install multi_key_dict

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
