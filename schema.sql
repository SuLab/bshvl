


drop schema if exists public cascade; create schema public;

DROP TABLE IF EXISTS genegene CASCADE;
CREATE TABLE genegene(
    docid text,
    mid1 text,
    mid2 text,
    mention1 text,
    mention2 text,
    is_correct boolean,
    features text[],
    sentence text,
    id bigint);

DROP TABLE IF EXISTS docids CASCADE;
CREATE TABLE docids(
    id bigint,
    docid text,
    folder text);

DROP TABLE IF EXISTS documents CASCADE;
CREATE TABLE documents(
    id bigint,
    docid text,
    document text);

DROP TABLE IF EXISTS sentences CASCADE;
CREATE TABLE sentences(
    id bigint,
    docid text,
    sentid text,
    sentence text,
    text text);

